.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(GenomicRanges)
library(parallel)
library(stringr)

nThread = 32
chrs = paste0("chr", c(1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
dir.create("data/MutTables/SomamutDB", showWarnings = F)


# load and filter sample Information table  #####
sampleInfo = read.table("data/rawdata/SomaMutDB_hg19/download.txt", sep = "\t", header = T, comment.char = "")
rownames(sampleInfo) = sampleInfo$Sample

# remove spaces from tissue names
sampleInfo$Tissue = gsub(x = sampleInfo$Tissue, pattern = " ", replacement = "_")
sampleInfo$Tissue..detail. = gsub(x = sampleInfo$Tissue..detail., pattern = " ", replacement = "_")


# for files where re-calling with SCcaller is available, remove the other file
# first, get all papers where there are Sccaller files available
SCcallerProjects = unique(sampleInfo$Paper[grep("S[c,C]caller", sampleInfo$Information)])
# go through them and manually remove non-SCcaller files
toRM = unlist(sapply(SCcallerProjects, function(proj){
  # project four: only remove MDA, keep the rest
  if(proj == "Aging and neurodegeneration are associated with increased mutations in single human neurons"){
    grep("MDA", sampleInfo$Sample[sampleInfo$Paper == proj], value = T, invert = T)
  }
  if(proj == "Single-cell analysis of somatic mutations in human bronchial epithelial cells in relation to aging and smoking"){
    c(grep("SCcaller", sampleInfo$Sample[sampleInfo$Paper == proj], value = T, invert = T),
      grep("PaperBulk", sampleInfo$Sample[sampleInfo$Paper == proj], value = T))
  }
  grep("SCcaller", sampleInfo$Sample[sampleInfo$Paper == proj], value = T, invert = T)
}))
sampleInfo = sampleInfo[!sampleInfo$Sample %in% toRM,]


# remove fetal samples/remove embryos (age = "Embryo")
sampleInfo = sampleInfo[sampleInfo$Age != "Embryo" | is.na(sampleInfo$Age),]


# remove samples marked as duplicates
# From https://www.science.org/doi/10.1126/science.aba8347#supplementary-materials:
# In order to test the reproducibility of variant calling from independent libraries for whole-
# genomes, two sets of duplicates (T08_61F_b01_lo0008 and T08_61F_b01_lo0048;
# T08_61F_b01_lo0035 and T08_61F_b01_lo0181) and one triplicate set
# (T08_61F_b01_lo0071, T08_61F_b01_lo0079 and T08_61F_b01_lo0091) were sequenced
# from microbiopsies of corresponding stretches of urothelium isolated from neighboring
# histology sections.
# --> Take file with most mutations each
toRM = c("T08_61F_b01_lo0048", "T08_61F_b01_lo0035",  "T08_61F_b01_lo0071", "T08_61F_b01_lo0091")
sampleInfo = sampleInfo[!sampleInfo$Sample %in% toRM,]

# I decided against removing some sutdies although they contain cancer patients etc. I will
# have to do a careful evaluation of the data (e.g. per study/patient) anyway.
#####


# split off exonic files and make sure they are not duplicated in WGS data #####
# WES is in sampleInfo$Information. Grep for WES and exonic
# WGS samples can be used for WES, but not vise versa
WESProjects = unique(sampleInfo$Paper[grep("WES|exonic", sampleInfo$Information)]) # 2554 samples are WES
# manually review contents:
# View(sampleInfo[sampleInfo$Paper == WESProjects[1],])
# 1, 2,3,5,6: all WES --> use for WES only
# 4: some sample WES and WGS, split accordingly using "WES" tag in Information column
# tag each sample in WES only, WGS only (when duplicated sample with WGS), or both
# create variable indicating whether we can use each sample for WGS or WES
sampleInfo$WES = T
sampleInfo$WGS = T
for(proj in WESProjects){
  subDat = sampleInfo[sampleInfo$Paper == proj,]
  subSamps = subDat$Sample[grep("WES|exonic",subDat$Information)]
  sampleInfo[subSamps,"WGS"] = F # mark exonic samples
  subSamps = subDat$Sample[grep("WES|exonic",subDat$Information, invert = T)]
  sampleInfo[subSamps,"WES"] = F # mark genomic samples
}

# Paper: Effects of psoriasis and psoralen exposure on the somatic mutation landscape of the skin:
# all WES (except: whole-genome sequencing (WGS) of 16 microbiopsies from three patients (patients 18,21 and 34))
# could not find info on which samples/muts affected, so use for WEX only, to be sure)
# --> mark all as WES only
sampleInfo$WGS[sampleInfo$Paper == "Effects of psoriasis and psoralen exposure on the somatic mutation landscape of the skin"] = F
# Substantial somatic genomic variation and selection for BCOR mutations in human induced pluripotent stem cells
# contains WES samples. indicated with _WES in sample ID
sampleIDsToProcess = sampleInfo$Sample[sampleInfo$Paper == "Substantial somatic genomic variation and selection for BCOR mutations in human induced pluripotent stem cells"]
seqType = data.frame(t(sapply(sampleIDsToProcess, function(x){
  if(length(grep("_WES|_hiWES", x))>0){
    c(x,strsplit(x, "_\\s*(?=[^_]+$)", perl=TRUE)[[1]])
  } else{
    c(x,x, "WGS")
  }
})))
seqTypeBySamp = split(seqType, seqType$X2)
# go through each sample: use only WGS sample for WGS. Use hiWES over WES if available.
for(samp in names(seqTypeBySamp)){
  tempDat = seqTypeBySamp[[samp]]
  if("WGS" %in% tempDat[,3]){ # remove all non-WGS files from WGS data
    sampleInfo[tempDat[tempDat$X3 != "WGS",1],"WGS"] = F
  }
  if("WES" %in% tempDat[,3]){
    sampleInfo[tempDat[tempDat$X3 == "WES",1],"WGS"] = F # remove this file from WGS, in case it wasn't done in the previous step yet
    sampleInfo[tempDat[tempDat$X3 != "WES",1],"WES"] = F # remove all other files from WES
  }
  if("hiWES" %in% tempDat[,3]){  # use hiWES over WES when available. don't use both
    sampleInfo[tempDat[tempDat$X3 == "hiWES",1],"WGS"] = F # remove this file from WGS, in case it wasn't done in the previous step yet
    sampleInfo[tempDat[tempDat$X3 != "hiWES",1],"WES"] = F # remove all other files from WES
  }
}

# The genomic landscapes of individual melanocytes from human skin
# most samples underwent WGS and WES, it is unclear which mutations come from where. Use for WES only, to be sure
# View(sampleInfo[sampleInfo$Paper == "The genomic landscapes of individual melanocytes from human skin",] )
sampleInfo$WGS[sampleInfo$Paper == "The genomic landscapes of individual melanocytes from human skin"] = F

#####

# create a Mega-table with all mutations plus annotation  #####
allMuts = sapply(sampleInfo$Sample, function(samp){
  # print(samp)
  if(samp == "T5_B11"){
    muts = read.table(paste0("data/rawdata/SomaMutDB_hg19/", samp, ".vcf"), quote = "", sep = "\t", fill = T)
  } else{
    muts = read.table(paste0("data/rawdata/SomaMutDB_hg19/", samp, ".vcf"), quote = "", sep = "\t")
  }
  muts = cbind(muts[,c(1,2,4,5)], samp)
  return(muts)
}, simplify = F)
dat = do.call(rbind,allMuts)
colnames(dat) = c("chr", "pos", "ref", "alt", "ID")
dat = cbind(dat, sampleInfo[dat$ID,])
rm(allMuts, sampleInfo)
gc()
#####



# filter #####
# only chrs 1-22
dat$chr = paste0("chr", dat$chr)
dat = dat[dat$chr %in% chrs,]


# Remove indels
indels = (nchar(dat$ref)!=1 | nchar(dat$ref)!=1)
dat = dat[!indels,]

#  filter for callable regions
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # gr

datGR = GRanges(seqnames=dat$chr,
                ranges=IRanges(start=dat$pos, end=dat$pos))
ov = is.na(findOverlaps(query=datGR, subject=gr, select = "first"))
dat = dat[ov,]
#####


# get sequence context for positions  #####
TPs_bed = data.frame(dat$chr, dat$pos-3, dat$pos+2)
TPs_bed = format(TPs_bed, trim = T)
tmpfile = tempfile(pattern = "healthyTissues", fileext = ".bed", tmpdir = "temp")
write.table(TPs_bed, file=tmpfile, quote=F, col.names=F, row.names=F, sep="\t")
cmd = paste0("bedtools getfasta -fi data/rawdata/GRCh37.primary_assembly.genome.fa ",
             "-bed  ",tmpfile,
             " -name -tab | cut -f 2")
context = system(cmd, intern = T)
dat$context = context
dat = dat[substring(context, 3,3) == dat$ref,]
file.remove(tmpfile)
save(dat, file = "data/MutTables/SomamutDB/SomamutDB_dat.RData")
#####

load("data/MutTables/SomamutDB/SomamutDB_dat.RData")
WEStissues = sort(unique(dat$Tissue[dat$WES]))
save(WEStissues, file = "data/MutTables/SomamutDB/WEStissues.RData")
WGStissues = sort(unique(dat$Tissue[dat$WGS]))
save(WGStissues, file = "data/MutTables/SomamutDB/WGStissues.RData")

# # WES data #####
# WESdat = dat[dat$WES,]
# # load sequences for TN generation
# exonSeqs = read.table("data/processedData/codExons_noUTR_filtered.withsequences.bed")
# load("data/processedData/codExons_noUTR_filtered.RData") # codExons_filtered
# codExons_filtered = GRanges(codExons_filtered)
# 
# # in this file, covered regions are extended by two bases in
# # both direction, in order to be able to find the 5mers.
# colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
# exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]
# # for each tissue
# set.seed(53634)
# dumpVar = sapply(sort(unique(WESdat$Tissue)), function(tissue){ 
#   print(tissue)
#   TPs = WESdat[WESdat$Tissue == tissue,]
# 
#   # store mutations to exclude during TN generation
#   posToFilter = paste(TPs$chr, TPs$pos, sep = "_")
# 
#   # exclude positions that were mutated more than once
#   dupl = duplicated(posToFilter)
#   TPs = TPs[!dupl,]
# 
# 
#   # subset to exome regions
#   TPsGR = GRanges(seqnames= TPs$chr,
#                   ranges=IRanges(start = TPs$pos,
#                                  end=TPs$pos))
#   ov = !is.na(findOverlaps(query=TPsGR, subject=codExons_filtered, select = "first"))
#   # print(table(ov))
#   TPs = TPs[ov,]
# 
# 
#   # Get TNs
#   TNs = lapply(chrs, function(chr){
#     # cat(chr,' ')
#     # subset data to chromosome
#     TPs_chr = TPs[TPs$chr == chr,]
#     exonSeqs_chr = exonSeqs[exonSeqs$chr == chr,]
#     pents = TPs_chr$context
#     TNs_chr = lapply(unique(pent2context[pents]), function(pent){
#       # cat(pent, ' ')
#       pent2search = fivemers[pent,]
#       nmatch1 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[1])
#       nmatch2 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[2])
#       nmatch = nmatch1 + nmatch2
#       weights = nmatch/sum(nmatch)
#       toRm = rep(1,sum(pents %in% pent2search))
#       pos = NULL
#       while(length(toRm) > 0){
#         TNs_sample = sample(1:nrow(exonSeqs_chr),
#                             size = length(toRm),
#                             prob = weights, replace = T)
#         new = t(sapply(TNs_sample, function(i){
#           # cat(i,' ')
#           gene = unlist(exonSeqs_chr[i,])
#           matches1 = cbind(gregexpr(pattern = pent2search[1], text = gene["sequence"])[[1]],pent2search[1])
#           matches2 = cbind(gregexpr(pattern = pent2search[2], text = gene["sequence"])[[1]],pent2search[2])
#           matches = rbind(matches1, matches2)
#           matches = matches[matches[,1] != "-1",, drop = F]
#           if(nrow(matches) == 1){
#             sampMatch = matches[1,]
#           } else{
#             sampMatch = matches[sample(1:nrow(matches),size = 1),]
#           }
#           samp = as.integer(sampMatch[1])
#           c(gene["chr"],
#             pos = (as.numeric(gene["start"]) + samp +2),
#             ref = substr(gene["sequence"],start = samp+2, stop = samp+2),
#             context = sampMatch[2])
#         }))
#         pos = rbind(pos,new)
#         temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
#         #  make sure that none of the positions are actually overlapping with existing positions:
#         toRm = which(temp %in% posToFilter | duplicated(temp))
#         if(length(toRm) > 0){
#           pos = pos[-toRm,]
#         }
#       }
#       return(pos)
#     })
#     TNs_chr = do.call(rbind,TNs_chr)
#   })
#   TNs = do.call(rbind,TNs)
#   TNs = as.data.frame(TNs)
#   TNs$pos = as.integer(TNs$pos)
#   colnames(TNs) = c("chr", "pos", "ref", "context")
# 
#   # save (also save meta information of TPs)
#   # save
#   TNs$alt = NA
#   TNs$mutated = 0
#   TPs$mutated = 1
#   meta = TPs
#   rownames(meta) = paste0(TPs$chr, "_", TPs$pos)
#   Muts = rbind(TPs[,c("chr", "pos", "ref", "alt", "context", "mutated")],
#                TNs[,c("chr", "pos", "ref", "alt", "context", "mutated")])
#   Muts = Muts[order(Muts[,1], Muts[,2]),]
#   ids = do.call(paste,c(Muts, sep = "_"))
#   MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]),
#                   ids)
#   save(Muts, meta, file = paste0("data/MutTables/SomamutDB/", tissue, "_WES.RData"))
# 
#   ids = do.call(paste,c(Muts, sep = "_"))
#   MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]),
#                   ids)
#   write.table(MutsBed, file = paste0("data/MutTables/SomamutDB/", tissue, "_WES.bed"),
#               col.names = F, row.names = F, sep = "\t", quote = F)
#   return(NULL)
# })
# 
# #####




# WGS data #####
load("data/MutTables/SomamutDB/SomamutDB_dat.RData")
load("data/MutTables/SomamutDB/WGStissues.RData")
WGSdat = dat[dat$WGS,]

# for each tissue
set.seed(24554)
dumpVar = sapply(sort(unique(WGSdat$Tissue)), function(tissue){
  print(tissue)
  TPs = WGSdat[WGSdat$Tissue == tissue,]
  
  # store mutations to exclude during TN generation
  posToFilter = paste(TPs$chr, TPs$pos, sep = "_")
  
  # exclude positions that were mutated more than once
  dupl = duplicated(posToFilter) 
  TPs = TPs[!dupl,]
  
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
      samp = sample(toSamp, size=n, replace = F)
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
  print("done with TNs")
  TNs = do.call(rbind,TNs)

  # save (also save meta information of TPs)
  TNs$ref = substr(TNs$context, 3,3)
  TNs$alt = NA
  TNs$mutated = 0
  TPs$mutated = 1
  
  meta = TPs
  rownames(meta) = paste0(TPs$chr, "_", TPs$pos)
  Muts = rbind(TPs[,c("chr", "pos", "ref", "alt", "context", "mutated")],
               TNs[,c("chr", "pos", "ref", "alt", "context", "mutated")])
  Muts = Muts[order(Muts[,1], Muts[,2]),]
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]),
                  ids)
  save(Muts,meta, file = paste0("data/MutTables/SomamutDB/", tissue, "_WGS.RData"))
  
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                  ids)
  write.table(MutsBed, file = paste0("data/MutTables/SomamutDB/", tissue, "_WGS.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)
  gc() 
  return(NULL)
})
#####
