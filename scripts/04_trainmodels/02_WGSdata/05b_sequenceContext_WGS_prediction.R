# preparation #####
dir.create("data/Modeling/WholeGenomeData/", showWarnings = F)
dir.create("data/Modeling/WholeGenomeData/contextModel", showWarnings = F)
tissues =  c("brain", "breast", "esophagus", "kidney", "liver", "ovary" ,
             "prostate", "skin", "allTissues")
load("data/Modeling/WholeGenomeData/contextModel/context_pPerChrWGS.RData")
load("data/processedData/pentamers.RData")
#####


# iterate through WGS parts and predict based on sequence context ####
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  pChr = pPerChrWGS[[tissue]]
  # use Final model for prediction
  p = pChr[pChr$chr == "all",]
  p = setNames(object=p$chrProb, nm=p$fivemers)
  
  sapply(1:22, function(i){ #iterate through chromosomes
    cat("chr",i, ' ')
    # load mutations
    load(paste0("data/MutTables/WholeGenomeData/Muts_WithContext_chr",i,".RData")) # Muts
    # position 10001 until 249240621
    cont = pent2context[Muts$context]
    predictions = data.frame(Muts$chr, Muts$pos, Muts$context, prediction = p[cont])
    save(predictions, 
         file = paste0("data/Modeling/WholeGenomeData/contextModel/WGSmuts_contextPredictions_chr",
                       i, "_", tissue, ".RData"))
    return(NA)
  })
  cat('\n')
  return(NA)
})
#####

