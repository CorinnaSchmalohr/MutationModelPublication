library(GenomicRanges)
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
# only once: go through all possible pentamers. use emboss to get pentamer positions #####
dir.create("data/processedData/pentLocations/", showWarnings=F)
# TODO filter for callable regions
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData")
library(parallel)
cl <- makeCluster(8, type="FORK")
start_time <- Sys.time()
dumpVar = parSapply(cl=cl, X=names(pent2context), FUN=function(pent){
  cat(pent, ' ')
  cmd = paste0("lib/EMBOSS-6.6.0/emboss/fuzznuc ",
               "-sequence data/rawdata/genome/GRCh37.primary_assembly.genome.fa ",
               "-pattern ", pent, " -rformat2 gff ",
               "-outfile temp/", pent, ".gff")
  system(cmd)
  temp = read.table(paste0("temp/", pent, ".gff"))
  pentLoc = GRanges(seqnames=temp$V1, 
                    ranges=IRanges(start=temp$V4, end=temp$V5))
  # pentLoc = setdiff(x = pentLoc, y = gr)
  pentLoc = subtract(x = pentLoc, y = gr, ignore.strand = T)
  pentLoc = unlist(pentLoc)
  pentLoc = pentLoc[width(pentLoc) == 5]
  save(pentLoc, file = paste0("data/processedData/pentLocations/", pent, ".RData"))
  file.remove(paste0("temp/", pent, ".gff"))
})
stopCluster(cl)
#####
finished = NA
save(finished, file = "temp/finished.RData")
