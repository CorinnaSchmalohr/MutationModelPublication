
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
library(readxl)
library(plyr)
nThreads = 28
tissues = c("brain", "breast", "esophagus", "kidney", "liver",
            "ovary", "prostate", "skin")

# get complete set of predictors #####
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000,
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
tab = read_xlsx("data/rawdata/dataMappingAlltissues_WGS.xlsx",
                sheet="allTissues", col_names=T)
tab$NA. = NULL
tab[tab == "NA"] = NA
tab = tab[,c(colnames(tab)[1:9],"allTissues")]
tab = tab[!is.na(tab[,"allTissues"]),]
# for predictors where we want multiple ranges, expand table
tab = apply(tab,1,function(x){
  if(is.na(x["range"])){
    return(x)
  } else{
    rangeWindow = strsplit(x["range"],";")[[1]]
    rangeInds = which(names(ranges) == rangeWindow[2]):which(names(ranges) == rangeWindow[1])
    subRanges = ranges[rangeInds]
    t(sapply(names(subRanges), function(r,y){
      y["range"] = subRanges[r]
      if(length(subRanges)>1){
        y["Name"] = paste0(y["Name"]," ",r)
        y["abbreviation"] = paste0(y["abbreviation"],"_",r)
      }
      return(y)
    },y=x))
  }
})
tab = do.call(rbind, tab)
tab = as.data.frame(tab)
predictors = tab$abbreviation
# # #####

# combine data from all tissues #####
# combine the data tables
combDat = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue,
                           "_mapped.RData")) # data
  dat = data.frame(data$pred,
                   mutated = as.factor(data$muts$mutated),
                   chr = data$muts$chr,
                   pos = data$muts$pos,
                   tissue = tissue)
  if(nrow(dat)>1e6){
    dat = dat[sample(1:nrow(dat), 1e6),] # cap the number of allowed mutations due to memory issues
  }
  return(dat)
}, simplify = F)
combDat = do.call(rbind.fill, combDat)
combDat = combDat[ ,colSums(is.na(combDat)) == 0] # removes all cols with NAs

# remove duplicates
positions = paste(combDat$chr, combDat$pos, sep = "_")
multiple = duplicated(positions) | duplicated(positions, fromLast=T)
combDat = combDat[!multiple,]
# combine the chromosome annotations
combChroms = combDat$chr

combTissue = combDat$tissue
combDat = combDat[,c(predictors[predictors %in% colnames(combDat)], "mutated")]
save(combDat, combChroms, combTissue,
     file = "data/MutTables/WholeGenomeData/WGSMuts_generalModel_Muts_mapped_processed.RData")
######

load("data/MutTables/WholeGenomeData/WGSMuts_generalModel_Muts_mapped_processed.RData")
chroms = unique(combChroms)


# grow forest with impurity_corrected #####
print("growing forests with impurity_corrected")
imp = sapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = combDat[combChroms != cr,]
  rf = ranger(mutated ~ ., data = trainData, importance = "impurity_corrected",
              write.forest = F, seed = 1234, num.threads =  nThreads,
              probability = T, verbose=F)
  chromimp = rf$variable.importance
  save(chromimp, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination",
                          "_importances_gini_",cr,".RData"))
  return(rf$variable.importance)
})
save(imp, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination",
                        "_importances_gini.RData"))
#####

# create rf for prediction and test on held-out chromosomes #####
print("predictions training")
dumpVar = lapply(chroms, function(cr){
  if(file.exists(paste0("data/Modeling/WholeGenomeData/RF/TissueCombination", "_",
                        cr, "_forPrediction.RData"))){
    cat(cr, " already done ")
    return(NA)
  }
  cat(cr, ' ')
  trainData = combDat[combChroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              probability = T, verbose=F)
  save(rf, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination", "_",
                         cr, "_forPrediction.RData"))
  rm(rf, trainData); gc()
  return(NA)
}); cat('\n')
print("predictions application")
predictions = lapply(chroms, function(cr){
  cat(cr, " ")
  load(paste0("data/Modeling/WholeGenomeData/RF/TissueCombination", "_",
                         cr, "_forPrediction.RData"))
  testData = combDat[combChroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  save(temp, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination", "_",
                           cr, "_prediction.RData"))
  return(temp)
}); cat('\n')

names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination",
                                "_predictions.RData"))
#####


# create final model #####
print("creating final model")
rf = ranger(mutated ~ ., data = combDat, importance = "impurity_corrected",
            write.forest = F, seed = 1234, num.threads =  nThreads,
            probability = T, verbose=F)
importance = rf$variable.importance
rm(rf); gc()
rf = ranger(mutated ~ ., data = combDat,
            write.forest = T, seed = 1234, num.threads =  nThreads,
            probability = T, verbose=F)
save(rf, importance, file = paste0("data/Modeling/WholeGenomeData/RF/TissueCombination_", 
                                   "finalModel.RData"))

#####

print("done")
