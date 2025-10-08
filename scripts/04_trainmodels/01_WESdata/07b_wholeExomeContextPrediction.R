
# preparation #####
dir.create("data/Modeling/WholeExomeData", showWarnings = F)
dir.create("data/Modeling/WholeExomeData/contextModel", showWarnings = F)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin", "allTissues")
load("data/Modeling/exomeTrainData/contextModel/context_pPerChr.RData")
load("data/processedData/pentamers.RData")
#####



# predict mutation based on sequence context for whole exome #####
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  pChr = pPerChr[[tissue]]
  sapply(1:50, function(i){
    cat(i, ' ')
    # load data
    load(paste0("data/MutTables/WholeExomeData/exomeMuts_part", i, ".RData")) # subMuts
    datpos = data.frame(chr = subMuts$chr, context = subMuts$context)
    chrs = unique(subMuts$chr)
    
    # Final model prediction
    p = pChr[pChr$chr == "all",]
    p = setNames(object=p$chrProb, nm=p$fivemers)
    cont = pent2context[datpos$context]
    predictions = data.frame(datpos, prediction = p[cont])

    save(predictions, 
         file = paste0("data/Modeling/WholeExomeData/contextModel/exomeMuts_contextPredictions_part",
                       i, "_", tissue, ".RData"))
    return(NA)
  })
  cat('\n')
  return(NA)
})

#####
