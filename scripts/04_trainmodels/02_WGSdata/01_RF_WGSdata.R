args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain", "breast", "esophagus", "kidney", "liver",
            "ovary", "prostate", "skin") 
tissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
nThreads = 28

dir.create("data/Modeling/WholeGenomeData/",showWarnings=F)
dir.create("data/Modeling/WholeGenomeData/RF",showWarnings=F)

print(tissue)

# load data for this tissue######
load(paste0("data/MutTables/WholeGenomeData/WGSMuts_", 
            tissue, "_mapped_processed.RData"))
chroms = unique(datchroms)
#####

# grow forest with impurity_corrected #####
print("growing forests with impurity_corrected")
imp = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = F, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F)
  res = rf$variable.importance
  save(res, file =  paste0("data/Modeling/WholeGenomeData/RF/", tissue,
                           "_importances_gini_", cr, ".RData"))
  rm(rf, trainData); gc()
  return(res)
})
print("saving")
save(imp, file = paste0("data/Modeling/WholeGenomeData/RF/", tissue,
                        "_importances_gini.RData"))
######


# predict on held-out chromosomes #####
print("prediction forest")
dumpVar = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = dat[datchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F)
  save(rf, file = paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_", 
                         cr, "_forPrediction.RData"))
  rm(rf, trainData); gc()
  return(NA)
})
print("predictions")
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_", 
                         cr, "_forPrediction.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  rm(rf, testData, p); gc()
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/WholeGenomeData/RF/", tissue,
                                "_predictions.RData"))
cat('\n')
#####



# create final model #####
rf = ranger(mutated ~ ., data = dat, importance = "impurity_corrected",
            write.forest = F, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
importance = rf$variable.importance
rf = ranger(mutated ~ ., data = dat,
            write.forest = T, seed = 1234, num.threads =  nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=F)
save(rf, importance, file = paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_", 
                                   "finalModel.RData"))

#####

print("done")
