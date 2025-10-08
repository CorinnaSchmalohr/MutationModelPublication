

# iterate through tissues and exome parts: load RF and context-based predictions
tissues = c( "brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin") # "allTissues"
capitalize <- function(x) {
  s <- strsplit(x, "_")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2),
        sep = "", collapse = " ")
}


dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  data = sapply(1:50, function(i){
    cat(i, " ")
    # load mutations
    load(paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,"_",
                tissue, "_mapped.RData")) # data
    # load context predictions: 
    load(paste0("data/Modeling/WholeExomeData/contextModel/exomeMuts_contextPredictions_part",
                i, "_", tissue, ".RData")) # predictions
    if(!all(predictions$context == data$muts$context)){
      warning(paste0("non-matching sequence context in tissue ", tissue, " part ", i))
    }
    res = cbind(data$muts, contextPreds = predictions$prediction)
    # load RF mutations
    load(paste0("data/Modeling/WholeExomeData/RF/exomeMuts_part",
                i, "_",tissue, "_RFpredictions.RData")) # predictions
    if(!all(predictions$context == data$muts$context)){
      warning(paste0("non-matching sequence context in tissue ", tissue, " part ", i))
    }
    res = cbind(res, RFpreds = predictions$prediction)
    
  }, simplify = F); cat("\n")
  data = do.call(rbind, data)
  data$mutated = NULL
  data$alt = NULL
  
  dataPos = paste(data$chr, data$pos, sep = "_")
  
  # annotate with cancer mutations 
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts.RData"))# Muts
  tempPos = paste(Muts$chr, Muts$pos, sep = "_")
  data$cancerMuts = as.integer(dataPos %in% tempPos)
  # annotate with healthy tissue mutations
  
  if(file.exists(paste0("data/MutTables/SomamutDB/", capitalize(tissue), "_WES.RData"))){
    load(paste0("data/MutTables/SomamutDB/", capitalize(tissue), "_WES.RData")) #Muts, meta
    tempPos = paste(Muts$chr, Muts$pos, sep = "_")
    data$healthyMuts = as.integer(dataPos %in% tempPos)
  } else
    data$healthyMuts = NA
  
  
  rm(meta, dataPos, Muts, tempPos); gc()
  
  # create simple combinations: multiply, mean, combOdds
  data$mult = data$RFpreds * data$contextPreds
  data$mean = (data$RFpreds + data$contextPreds)/2
  odds = data$RFpreds/(1-data$RFpreds) * data$contextPreds/(1-data$contextPreds)
  data$combOdds = odds/(odds+1)
  # linear Model: logistic regression taking RF preds and context preds as input and TP/TN as output.
  # with CWCV
  LMcombination = sapply(unique(data$chr), function(chrom){
    cat(chrom, ' ')
    trainData = data[data$chr != chrom,c("contextPreds", "RFpreds", "cancerMuts")]
    model = glm(formula = cancerMuts ~ RFpreds + contextPreds, data = trainData, 
                       family = binomial(link = "logit"), x=F, y =F, model=F)
    testData = data[data$chr == chrom,c("chr", "pos", "contextPreds", "RFpreds")]
    yhat = predict(model, newdata = testData, type = "response")
    return(yhat)
  }, simplify = F); cat("\n")
  data$LMcombination = do.call(c, LMcombination)
  save(data, file = paste0("data/Modeling/WholeExomeData/combinedPredictions/combinedPredictions_", tissue, ".RData"))
  return(NA)
})



