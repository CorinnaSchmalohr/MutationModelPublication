dir.create("data/Modeling/exomeTrainData/contextModel/", showWarnings = F)
# preparation #####
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin", "allTissues")
library(Biostrings)
load("data/processedData/pentamers.RData")



# get baseline pentamer frequencies per chromosome #####
load("data/MutTables/WholeExomeData/exomeMuts.RData")
chrFreqs = sapply(c(paste0("chr", 1:22), "all"), function(cr){
   if(cr == "all"){
     cat("Whole data")
     sub = Muts$context
     fr = table(factor(sub, levels = fivemers[,1]))
     fr2 = table(factor(sub, levels = fivemers[,2]))
     fr+fr2
   } else {
    cat(cr, ' ')
    sub = Muts$context[Muts$chr == cr]
    fr = table(factor(sub, levels = fivemers[,1]))
    fr2 = table(factor(sub, levels = fivemers[,2]))
    fr+fr2
   }
   
}) 
save(chrFreqs, file = "data/Modeling/exomeTrainData/contextModel/chrFreqs.RData")
#####
# load("data/rdata/chrFreqs.RData")



# get probability of mutation based on sequence context #####
pPerChr = sapply(tissues, function(tissue){
   print(tissue)
   # load(paste0("data/rdata_deprecated/", tissue, "/Muts.RData"))
   load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts.RData"))

   temp = sapply(c(paste0("chr", 1:22), "all"), function(cr){
      if(cr == "all"){
        # Context model on whole data ?
        sub = Muts$context[Muts$mutated == 1]
      } else {
        sub = Muts$context[Muts$chr != cr & Muts$mutated == 1]
      }
      fr1 = table(factor(sub, levels=fivemers[,1]))
      fr2 = table(factor(sub, levels=fivemers[,2]))
      fr = as.numeric(fr1+fr2)
      res = data.frame(chr = cr,fivemers = names(fr1),obs = fr,
                      chrFreq = rowSums(chrFreqs) - chrFreqs[,cr],
                      chrProb = fr/(rowSums(chrFreqs) - chrFreqs[,cr]))
      return(res)
   }, simplify=F)
   temp = do.call(rbind,temp)
   temp$tissue = tissue
   return(temp)
}, simplify = F)

save(pPerChr, file = "data/Modeling/exomeTrainData/contextModel/context_pPerChr.RData")
######


