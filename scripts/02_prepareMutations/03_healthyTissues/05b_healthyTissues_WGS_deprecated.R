args = as.numeric(commandArgs(trailingOnly=T))
tissues =c("adipocytes", "adrenal_gland", "bladder", 
           "bonemarrow", "brain", "breast", "colon", 
            "esophagus", "fibroblast", "heart",  "kidney",
           "liver", "lung", "pancreas", "placenta", "prostate", "rectum", 
           "skeletal_muscle", "skin", "small_intestine", "spleen", 
           "stomach", "testicle", "thyroid", "tonsil", "ureter") #"iPSC","blood", "embryonic_stem_cell", "endometrium",
tissue = tissues[args]
print(tissue)


# Packages ######
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(GenomicRanges)
library(parallel)
######

# preparation #####
nThread = 10
dir.create("data/MutTables/healthyTissuesWGS/", showWarnings = F)
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # gr


chrs = paste0("chr", c(1:22))
chrTranslator = setNames(c(chrs, chrs),
                         c(chrs, 1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
#####

#  prepare data #####

TPs = read.table(paste0("data/rawdata/SomaMutDB/",tissue,".tsv"),
                 sep = "\t", header = T, comment.char = "")

# Remove indels 
TPs <- TPs[nchar(TPs$ref) == 1 & nchar(TPs$alt) == 1,]

# only chrs 1-22
colnames(TPs)[colnames(TPs) == "X.chr"] = "chr"
TPs$chr = chrTranslator[TPs$chr]
TPs = TPs[!is.na(TPs$chr),]

# store mutations to exclude during TN generation
posToFilter = paste(TPs$chr, TPs$pos, sep = "_")

# exclude positions that were mutated more than once
dupl = duplicated(posToFilter) 
TPs = TPs[!dupl,]

# get sequence context for positions
TPs_bed = data.frame(TPs$chr, TPs$pos-3, TPs$pos+2)
TPs_bed = format(TPs_bed, trim = T)
tmpfile = tempfile(pattern = tissue, fileext = ".bed", tmpdir = "temp")
write.table(TPs_bed, file=tmpfile, quote=F, col.names=F, row.names=F, sep="\t")
cmd = paste0("bedtools getfasta -fi data/rawdata/GRCh37.primary_assembly.genome.fa ", 
             "-bed  ",tmpfile,
             " -name -tab | cut -f 2")
context = system(cmd, intern = T)
TPs$context = context
if(!all(substring(context, 3,3) == TPs$ref)){
  TPs = TPs[substring(context, 3,3) == TPs$ref,]
}
file.remove(tmpfile)

#  filter for callable regions
TPsGR = GRanges(seqnames=TPs$chr, 
                ranges=IRanges(start=TPs$pos, end=TPs$pos))
ov = is.na(findOverlaps(query=TPsGR, subject=gr, select = "first"))
TPs = TPs[ov,c("chr", "pos", "ref", "alt","context")]

# Get TNs
cl <- makeCluster(nThread, type="FORK")
TNs = parLapply(cl = cl, unique(pent2context[TPs$context]), function(pent){
  pent2search = fivemers[pent,]
  samplePos = do.call(c,lapply(pent2search, function(p){
    load(paste0("data/processedData/pentLocations/", p, ".RData"))
    pentLoc$context = p
    return(pentLoc)
  }))
  start(samplePos) = start(samplePos) + 2
  width(samplePos) = 1
  toFilter = paste(seqnames(samplePos), start(samplePos), sep = "_")
  samplePos = samplePos[!toFilter %in% posToFilter]
  TNchr = lapply(as.character(unique(TPs$chr)), function(cr){ 
    n = sum(TPs$context[TPs$chr == cr] %in% pent2search)
    toSamp = which(seqnames(samplePos) == cr)
    samp = sample(toSamp, size=n)
    return(samplePos[samp])
  })
  TNchr = do.call(c, TNchr)
  TNchr = as.data.frame(TNchr)
  TNchr = TNchr[,c("seqnames", "start", "context")]
  colnames(TNchr) = c("chr", "pos", "context")
  TNchr$chr = as.character(TNchr$chr)
  return(TNchr)
})
stopCluster(cl)
TNs = do.call(rbind,TNs)

# save
TNs$ref = substr(TNs$context, 3,3)
TNs$alt = NA
TNs$mutated = 0
TPs$mutated = 1
Muts = rbind(TNs,TPs)
Muts = Muts[order(Muts[,1], Muts[,2]),]
Muts = Muts[,c("chr", "pos", "ref", "alt", "context", "mutated")]
ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
save(Muts, file = paste0("data/MutTables/healthyTissuesWGS/", tissue, "_SomaMutDB.RData"))

ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
write.table(MutsBed, file = paste0("data/MutTables/healthyTissuesWGS/", tissue, "_SomaMutDB.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
#####