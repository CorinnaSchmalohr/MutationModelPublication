# preparation
args = as.numeric(commandArgs(trailingOnly=T))
tissue2Cancer = c("brain" =  "CNS-Medullo\\|CNS-PiloAstro",
                  "breast" =  "Breast-AdenoCa",
                  "esophagus" =  "Eso-AdenoCa",
                  "kidney" =  "Kidney-RCC",
                  "liver" =  "Liver-HCC",
                  "ovary" =  "Ovary-AdenoCA",
                  "prostate" = "Prost-AdenoCA",
                  "skin" = "Skin-Melanoma") # colon and luad missing
tissue = names(tissue2Cancer)[args]
print(tissue)
nThread =16
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(GenomicRanges)
library(parallel)

chrs = paste0("chr", c(1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
#####

# get the column names
cmd = paste0("zcat -f data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz 2>/dev/null | head -n1 - ")
colnames = colnames(read.table(text = system(cmd, intern = T, ignore.stderr=), header = T, sep = "\t", quote=""))

# extract mutations for this tissue
cmd = paste0("zcat data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz | grep '",tissue2Cancer[tissue], "'")
muts = read.table(text = system(cmd, intern = T), header = F, sep = "\t", quote="")
colnames(muts) = colnames

print("subset to autosomic chromosomes")
muts$Chromosome = paste0("chr", muts$Chromosome)
table(muts$Chromosome %in% chrs, useNA = "always")
muts = muts[muts$Chromosome %in% chrs,]

print("remove indels")
table(muts$Variant_Type)
muts = muts[muts$Variant_Type == "SNP",]

# extract positions to exclude during TN sampling
posToFilter = cbind(muts$Chromosome,
                    paste(muts$Chromosome, muts$Start_position, sep = "_"))
save(posToFilter, file = paste0("data/MutTables/WholeGenomeData/posToFilter_", 
                                tissue, ".RData"))

print("filter out variants with too little depth/poor quality")
table(!is.na(muts$t_ref_count) &
        !is.na(muts$t_alt_count) &
        muts$t_alt_count >= 2 &
        muts$t_alt_count+muts$t_ref_count >=14 &
        muts$i_NumCallers >= 2, useNA = "always")
muts = muts[!is.na(muts$t_ref_count) & 
              !is.na(muts$t_alt_count) & 
              muts$t_alt_count >= 2 & 
              muts$t_alt_count+muts$t_ref_count >=14 & 
              muts$i_NumCallers >= 2,]

print("remove variants with >1% population frequency")
temp = muts$i_1000genomes_AF
temp[temp == ""] = NA
temp[grep("|", temp, fixed = T)] =  sapply(strsplit(temp[grep("|", temp, fixed = T)], "|", fixed = T), 
                                           function(x){max(as.numeric(x), na.rm = T)})
temp = as.numeric(temp)
table(is.na(temp) | temp < 0.01, useNA = "always")
muts = muts[is.na(temp) | temp < 0.01,]


print("exclude positions that were mutated more than once")
positions = paste(muts$Chromosome, muts$Start_position, sep = "_")
multiple = duplicated(positions) | duplicated(positions, fromLast=T)
table(multiple)
muts = muts[!multiple,]



#  filter for callable regions
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # this makes sense since we integrate with ChiP-seq data
muts$pentContext = toupper(substr(muts$ref_context, 9,13))
TPsGR = GRanges(seqnames=muts$Chromosome, 
              ranges=IRanges(start=muts$Start_position, end=muts$End_position), 
              context = muts$pentContext)
ov = is.na(findOverlaps(query=TPsGR, subject=gr, select = "first"))
TPs = muts[ov,c("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2", "pentContext")]
colnames(TPs) = c("chr", "pos", "ref", "alt","context")
save(TPs, file = paste0("data/MutTables/WholeGenomeData/TPs_", tissue,".RData"))
counts = table(muts$Tumor_Sample_Barcode[ov])
save(counts, file = paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue, "_nPerSample.RData"))
if(nrow(TPs)>3.5e6){
  TPs = TPs[sample(1:nrow(TPs), 3.5e6),]
}

# get TNs #####
# sample from these positions. 
print("getting TNs")
cl <- makeCluster(nThread, type="FORK")
TNs = parLapply(cl = cl, unique(pent2context[TPs$context]), function(pent){
  cat(pent, ' ')
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
}); cat('\n')
print("done with TNs")
stopCluster(cl)
TNs = do.call(rbind,TNs)
#####

# save #####
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
save(Muts, file = paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue, ".RData"))
write.table(MutsBed, file = paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue, ".bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
print("end")
#####
# })


