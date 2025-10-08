.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(GenomicRanges)
library(parallel)
nThread = 16
tissues = c("brain", "breast", "esophagus", "kidney", "liver",
            "ovary", "prostate", "skin") 

# load TPs and combine ##### 
combTPs = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/WholeGenomeData/TPs_", tissue,".RData"))
  if(nrow(TPs)>2e6){
    TPs = TPs[sample(1:nrow(TPs), 2e6),] # cap the number of allowed mutations due to memory issues
  }
  return(TPs)
}, simplify = F)
# sum(sapply(combTPs, nrow))
combTPs = do.call(rbind,combTPs)
print("exclude positions that were mutated more than once")
positions = paste(combTPs$chr, combTPs$pos, sep = "_")
multiple = duplicated(positions) | duplicated(positions, fromLast=T)
table(multiple)
combTPs = combTPs[!multiple,]
#####

# get TNs #####
posToFilterTotal = do.call(rbind,sapply(tissues, function(tissue){
  load( paste0("data/MutTables/WholeGenomeData/posToFilter_", 
             tissue, ".RData")) #posToFilter
  return(posToFilter)
}))
chrs = paste0("chr", c(1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
# get TNs #####
# sample from these positions. 
print("getting TNs")
cl <- makeCluster(nThread, type="FORK")
TNs = parLapply(cl = cl, unique(pent2context[combTPs$context]), function(pent){
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
  samplePos = samplePos[!toFilter %in% posToFilterTotal[,2]]
  TNchr = lapply(as.character(unique(combTPs$chr)), function(cr){ 
    n = sum(combTPs$context[combTPs$chr == cr] %in% pent2search)
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
combTPs$mutated = 1
Muts = rbind(TNs,combTPs)
Muts = Muts[order(Muts[,1], Muts[,2]),]
Muts = Muts[,c("chr", "pos", "ref", "alt", "context", "mutated")]
ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
save(Muts, file = paste0("data/MutTables/WholeGenomeData/WGSMuts_allTissues.RData"))
write.table(MutsBed, file = paste0("data/MutTables/WholeGenomeData/WGSMuts_allTissues.bed"),
            col.names = F, row.names = F, sep = "\t", quote = F)
print("end")
#####

