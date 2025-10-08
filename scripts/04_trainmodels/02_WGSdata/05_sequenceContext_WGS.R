dir.create("data/Modeling/WholeGenomeData/contextModel/", showWarnings = F)
# preparation #####
tissues =  c("brain", "breast", "esophagus", "kidney", "liver", "ovary" ,
             "prostate", "skin", "allTissues")
library(Biostrings)
library(GenomicRanges)
load("data/processedData/pentamers.RData") # fivemers, pent2context
chrs = paste0("chr", 1:22)
#####


# get baseline pentamer frequencies per chromosome #####
chrFreqsWGS = apply(fivemers,1, function(pents){
  load(paste0("data/processedData/pentLocations/", pents[1], ".RData")) # pentLoc
  count1 = table(factor(seqnames(pentLoc), levels = chrs))
  load(paste0("data/processedData/pentLocations/", pents[2], ".RData")) # pentLoc
  count2 = table(factor(seqnames(pentLoc), levels = chrs))
  return(count1 + count2)
})
chrFreqsWGS = t(chrFreqsWGS)
chrFreqsWGS = cbind(chrFreqsWGS, "all" = rowSums(chrFreqsWGS))
save(chrFreqsWGS, file = "data/Modeling/WholeGenomeData/contextModel/chrFreqsWGS.RData")
#####



# get probability of mutation based on sequence context #####
pPerChrWGS = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue, ".RData"))
  temp = sapply(c(paste0("chr", 1:22), "all"), function(cr){
    if(cr == "all"){
      # Context model on whole data 
      sub = Muts$context[Muts$mutated == 1]
    } else {
      sub = Muts$context[Muts$chr != cr & Muts$mutated == 1]
    }
    fr1 = table(factor(sub, levels=fivemers[,1]))
    fr2 = table(factor(sub, levels=fivemers[,2]))
    fr = as.numeric(fr1+fr2)
    res = data.frame(chr = cr,fivemers = names(fr1),obs = fr,
                     chrFreq = rowSums(chrFreqsWGS) - chrFreqsWGS[,cr],
                     chrProb = fr/(rowSums(chrFreqsWGS) - chrFreqsWGS[,cr]))
    return(res)
  }, simplify=F)
  temp = do.call(rbind,temp)
  temp$tissue = tissue
  return(temp)
}, simplify = F)

save(pPerChrWGS, file = "data/Modeling/WholeGenomeData/contextModel/context_pPerChrWGS.RData")
######


