
library(GenomicRanges)
tissues =  c("brain", "breast", "esophagus", "kidney", "liver", "ovary" ,
             "prostate", "skin")
capitalize <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}
dir.create("data/Modeling/WholeGenomeData/combinedPredictions", showWarnings = F)


# iterate through tissues and genome chromosomes parts: load RF and context-based predictions
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  # load cancer mutations 
  load(paste0("data/MutTables/WholeGenomeData/TPs_", tissue,".RData")) # TPs
  # cancerPos = paste(TPs$chr, TPs$pos, sep = "_")
  cancerGR = GRanges(seqnames=TPs$chr, 
                     ranges=IRanges(start=TPs$pos, width = 1))
  rm(TPs)
  # load healthy mutations
  if(file.exists(paste0("data/MutTables/SomamutDB/", capitalize(tissue), "_WGS.RData"))){
    load(paste0("data/MutTables/SomamutDB/", capitalize(tissue), "_WGS.RData")) #Muts, meta
    healthyGR = GRanges(seqnames=Muts$chr, 
                        ranges=IRanges(start=Muts$pos, width = 1))
    rm(Muts, meta)
  } 
  gc()
  i = 22 # use only chromosome 1 for now
  # dumpVar2 = sapply(1:22, function(i){
    # cat(i, " ")
    cancerGR = cancerGR[seqnames(cancerGR) == paste0("chr",i)]
    if(tissue != "ovary"){
      healthyGR = healthyGR[seqnames(healthyGR) == paste0("chr",i)]
    }
    # load context predictions: 
    load(paste0("data/Modeling/WholeGenomeData/contextModel/WGSmuts_contextPredictions_chr",
                i, "_", tissue, ".RData")) # predictions
    # nrow = 225041628
    data = GRanges(seqnames=predictions$Muts.chr, 
                   ranges=IRanges(start=predictions$Muts.pos, width = 1))
    data$contextPreds = predictions$prediction
    rm(predictions); gc()
    
    # load RF predictions part 1
    load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_MutsResult_chr",
                i,"_parts_1_50_RFpredictions.RData")) # predictions
    RFpred1 = predictions
    # part 2
    load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_MutsResult_chr",
                i,"_parts_51_100_RFpredictions.RData")) # predictions
    # RFpred = rbind(RFpred1, predictions)
    # TODO There are duplicated positions in the RF predictions. 
    dupl = duplicated(c(RFpred1$pos, predictions$pos))
    data$RFpreds = c(RFpred1$prediction, predictions$prediction)[!dupl]
    rm(predictions, RFpred1); gc()
    data
    
    # annotate with cancer and healthy mutations 
    cancerMatch = findOverlaps(query = data, subject = cancerGR, select = "first")
    data$cancerMuts = as.integer(!is.na(cancerMatch))
    if(tissue != "ovary"){
      healthyMatch = findOverlaps(query = data, subject = healthyGR, select = "first")
      data$healthyMuts = as.integer(!is.na(healthyMatch))
    }
    
    # create simple combinations: multiply, combOdds
    data$mult = data$RFpreds * data$contextPreds
    odds = data$RFpreds/(1-data$RFpreds) * data$contextPreds/(1-data$contextPreds)
    data$combOdds = odds/(odds+1)
    
    save(data, file = paste0("data/Modeling/WholeGenomeData/combinedPredictions/combinedPredictions_", 
                             tissue, "_chr",i,".RData"))
    
    # return(NA)
  # }, simplify = F); cat("\n")
  return(NA)
})



