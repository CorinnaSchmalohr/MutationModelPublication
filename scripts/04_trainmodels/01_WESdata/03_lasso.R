args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
dir.create("data/Modeling/exomeTrainData/Lasso",showWarnings=F)
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(c060)
print(tissue)

# prepare data for this tissue######
load(paste0("data/MutTables/exomeTrainData/", 
            tissue, "_Muts_mapped_processed.RData"))
chroms = unique(datchroms)
#####

# CWCV stability selection #####
print("CWCV")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  
  print("training model")
  trainData = dat[datchroms != cr,]
  trainX = as.matrix(trainData[,colnames(trainData) != "mutated"])
  trainY = as.numeric(levels(trainData$mutated))[trainData$mutated]
  sp = stabpath(x = trainX, y = trainY, family = "binomial")
  stab = stabsel(sp, error = 0.1, type = "pcer", pi_thr = 0.6)
  
  save(sp, stab, file = paste0("data/Modeling/exomeTrainData/Lasso/",
                               tissue, "_", cr, ".RData"))
  
  print("training significant models")
  sigFeatures = names(stab$stable)
  sigTrain = cbind(trainData[,sigFeatures], mutated = trainData$mutated)
  logR = glm(formula = mutated ~ ., data = sigTrain, 
             family = binomial(link = "logit"), x=F, y =F, model=F)
  save(logR, sigFeatures, file = paste0("data/Modeling/exomeTrainData/Lasso/",
                                        tissue, "_", cr, "_sig.RData"))
  
  
  print("predictions on left out chromosome")
  testData = dat[datchroms == cr,c(sigFeatures, "mutated")]
  yhat = predict(logR, newdata = testData, type = "response")
  temp = data.frame(pred  = yhat,label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/exomeTrainData/Lasso/", 
                                tissue, "_predictions_sig.RData"))
cat('\n')

# List all CWCV significant predictors
feature_selection = sapply(chroms, function(cr){
  load(paste0("data/Modeling/exomeTrainData/Lasso/", tissue, "_", cr, "_sig.RData"))
  return(sigFeatures)
})
names(feature_selection) = chroms
save(feature_selection, file = paste0("data/Modeling/exomeTrainData/Lasso/", 
                                      tissue,"_CWCV_feature_selection.RData"))
#####


# Final glm #####
print("training final glm")
X = as.matrix(dat[,colnames(dat) != "mutated"])
Y = as.numeric(levels(dat$mutated))[dat$mutated]
sp = stabpath(x = X, y = Y, family = "binomial")
stab = stabsel(sp, error = 0.1, type = "pcer", pi_thr = 0.6)
save(sp, stab, file = paste0("data/Modeling/exomeTrainData/Lasso/", tissue, ".RData"))

sigFeatures = names(stab$stable)

print("training final model on significant features")
sigTrain = cbind(dat[,sigFeatures], mutated = dat$mutated)
logR = glm(formula = mutated ~ ., data = sigTrain, 
           family = binomial(link = "logit"), x=F, y =F, model=F)
save(logR, sigFeatures, file = paste0("data/Modeling/exomeTrainData/Lasso/",
                                      tissue, "_sig.RData"))
#####
