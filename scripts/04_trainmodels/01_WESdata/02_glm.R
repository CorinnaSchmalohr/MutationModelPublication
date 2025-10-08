# args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin") # 
# tissue = tissues[args]
dir.create("data/Modeling/exomeTrainData/GLM",showWarnings=F)



dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  
  # prepare data for this tissue######
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped_processed.RData"))
  chroms = unique(datchroms)
  #####
  
  # CWCV glm #####
  print("CWCV")
  predictions = lapply(chroms, function(cr){
    cat(cr, ' ')
    
    print("training model")
    trainData = dat[datchroms != cr,]
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"), x=F, y =F, model=F)
    save(logR, file = paste0("data/Modeling/exomeTrainData/GLM/",
                             tissue, "_", cr, ".RData"))
    
    print("p-values")
    pvals = coef(summary(logR))[,4][-1]
    sigFeatures = names(pvals[pvals < 0.05])
    
    print("training significant models")
    trainData = cbind(trainData[sigFeatures], mutated = trainData$mutated)
    logR = glm(formula = mutated ~ ., data = trainData, 
               family = binomial(link = "logit"), x=F, y =F, model=F)
    save(logR, sigFeatures, file = paste0("data/Modeling/exomeTrainData/GLM/", 
                                          tissue, "_", cr, "_sig.RData"))
    
    print("predictions on left out chromosome")
    testData = cbind(dat[datchroms == cr,sigFeatures], mutated = dat[datchroms == cr,]$mutated)
    yhat = predict(logR, newdata = testData, type = "response")
    temp = data.frame(pred  = yhat,label = testData$mutated)
    return(temp)
  })
  names(predictions) = chroms
  save(predictions, file = paste0("data/Modeling/exomeTrainData/GLM/", tissue, 
                                  "_predictions_sig.RData"))
  cat('\n')
  
  # List all CWCV significant predictors
  feature_selection = sapply(chroms, function(cr){
    cat(cr, ' ')
    load(paste0("data/Modeling/exomeTrainData/GLM/", tissue, "_", cr, "_sig.RData"))
    return(sigFeatures)
  })
  names(feature_selection) = chroms
  save(feature_selection, file = paste0("data/Modeling/exomeTrainData/GLM/", 
                                        tissue,"_CWCV_feature_selection.RData"))
  cat('\n')
  #####
  
  
  # Final glm #####
  print("training final glm")
  logR = glm(formula = mutated ~ ., data = dat, family = binomial(link = "logit"),
             x=F, y =F, model=F)
  save(logR, file = paste0("data/Modeling/exomeTrainData/GLM/", tissue, ".RData"))
  
  print("p-values")
  pvals = coef(summary(logR))[,4][-1]
  sigFeatures = names(pvals[pvals < 0.05])
  # check if sig. features are also present in at least 60% of the selected CWCV features
  # feature_selection = table(do.call(c,feature_selection))
  # feature_selection = feature_selection[feature_selection > length(chroms)*0.6]
  # sigFeatures_selection = sigFeatures[sigFeatures %in% names(feature_selection)]
  
  print("training final model on significant features")
  trainData = cbind(dat[sigFeatures], mutated = dat$mutated)
  logR = glm(formula = mutated ~ ., data = trainData, 
             family = binomial(link = "logit"), x=F, y =F, model=F)
  save(logR, sigFeatures, 
       file = paste0("data/Modeling/exomeTrainData/GLM/", tissue, "_sig.RData"))
  #####
  
  print("done with this tissue")
  
})

randomVar = "finished"
save(randomVar, file = "temp/finishedWithGLM.R")