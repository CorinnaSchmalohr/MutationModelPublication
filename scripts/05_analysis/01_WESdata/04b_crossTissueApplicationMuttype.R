.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)
library(ROCR)
library(ggplot2)
source("scripts/05_analysis/00_NamesAndColors.R")
cr= "chr1"
load("data/processedData/dataInfos.RData")
mindatsize = min(sapply(dataInfos, function(x)x$nMuts))
dir.create("data/Modeling/exomeTrainData/CrossTissue", showWarnings = F)
dir.create("fig/CrossTissue", showWarnings = F)
plotEnding = "_20230902"
colorRange = c(0.45,0.7) 
nThreads = 32
# prepare mutation type and fivemers #####
mutClass = rbind(c("T>A", "A>T"),
                 c("T>C", "A>G"),
                 c("T>G", "A>C"),
                 c("C>A", "G>T"),
                 c("C>G", "G>C"),
                 c("C>T", "G>A"))
rownames(mutClass) = apply(mutClass, 1,paste, collapse = "/")
mutClassTranslator = setNames(nm=c(mutClass[,1], mutClass[,2]),
                              object=c(rownames(mutClass), rownames(mutClass)))
load("data/processedData/pentamers.RData")
fivemerTranslator  = setNames(nm=c(fivemers[,1], fivemers[,2]),
                              object=c(rownames(fivemers), rownames(fivemers)))
#####


# stratify cross-tissue analysis by substitution type #####
# Take the full model (all mutation types) from each tissue and use it to predict 
# for the data from another tissue (I already have this)

performanceCrossTissueSubstType = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  #   load data for tissue2predict --> testData
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue2predict, "_Muts_mapped.RData")) # data
  testData = cbind(data$pred, mutated = as.factor(data$muts$mutated))
  testData = as.data.frame(testData)
  testAnnot = data.frame(chr = data$muts$chr,
                         mutated = testData$mutated,
                         context = data$muts$context,
                         substitution = mutClassTranslator[paste0(data$muts$ref, ">", data$muts$alt)])
  # infer mutationtype for TNs based on context
  for(context in rownames(fivemers)){
    testAnnot$substitution[(testAnnot$context %in% fivemers[context,]) & testAnnot$mutated == 0] = 
      sample(testAnnot$substitution[(testAnnot$context %in% fivemers[context,]) & testAnnot$mutated == 1])
  }
  crossTissuePerfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    # load data --> trainData
    load(paste0("data/MutTables/exomeTrainData/",
                predictingTissue, "_Muts_mapped_processed.RData")) # dat, datchroms
    trainData = dat
    trainChroms = datchroms
    rm(dat, datchroms)
    # match predictors between test and trainData
    predictors = colnames(trainData)[-(ncol(trainData))]
    testDataMatched = testData[,colnames(testData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testDataMatched)]
    for(x in toAdd){ # Add mean value of predictors missing for the tissue to predict on (from predicting tissue)
      testDataMatched[x] = rep(mean(trainData[,x]), nrow(testDataMatched))
    }
    # iterate through chromosomes (for CWCV)
    chromWisePreds = sapply(names(chrCols), function(tempCr){
      cat(tempCr, ' ')
      # load model from predictingTissue
      load(paste0("data/Modeling/exomeTrainData/RF/", predictingTissue, "_", 
                  tempCr, "_forPrediction.RData")) # rf
      # subset testData to chromosome 
      testDataChrom = testDataMatched[testAnnot$chr == tempCr,]
      # predict testData
      yHat = predict(rf,data=testDataChrom,type="response", num.threads=nThreads)
      # cbind predictions, trueMutationState, context, mutationtype --> predictionResult
      res = cbind(prediction = yHat$predictions[,2], 
                  testAnnot[testAnnot$chr == tempCr,])
      # return predictions
      return(res)
    }, simplify = F)
    # rbind predictions
    preds = do.call(rbind, chromWisePreds)
    # for subst in substitutions{
    perfBySubstType = sapply(rownames(mutClass), function(subst){
      #   subset predictionResult to fivemer
      subPreds = preds[preds$substitution == subst,]
      #   compute AUC
      perf = prediction(subPreds$prediction, subPreds$mutated)
      auc = performance(perf,"auc")@y.values[[1]]
      #   return AUC
      return(auc)
    })
    resTable = data.frame(trainTissue = predictingTissue, 
                     testTissue = tissue2predict,
                     substitution = names(perfBySubstType),
                     AUC = perfBySubstType, row.names = NULL)
    return(resTable)
  }, simplify = F); cat("\n")
  crossTissuePerfs = do.call(rbind,crossTissuePerfs)
  return(crossTissuePerfs)
}, simplify = F)    
save(performanceCrossTissueSubstType, 
     file = "data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissueSubstType.RData")

# Visualize using side-by-side heatmaps.
# load("data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissueSubstType.RData")
# plotDat = do.call(rbind,performanceCrossTissueSubstType)
# plotDat$trainTissue = paste0(t2T[plotDat$trainTissue], " Model")
# plotDat$testTissue = t2T[plotDat$testTissue]
# plotDat$substitution = mutClass[plotDat$substitution,1]
# ggplot(plotDat, aes(substitution,testTissue , fill= AUC)) +
#   geom_tile() +
#   # scale_fill_viridis_c(limits=colorRange) +
#   scale_fill_viridis_c()+
#   theme_minimal() +
#   facet_wrap(~trainTissue, nrow = 2, scales = "free_x") +
#   labs(x = "Substitution", y = "Test data") 
# ggsave(filename = paste0("fig/CrossTissue/crossTissueSubstitutionType", 
#                          plotEnding, ".png"),
#        width = 25, height = 15, units = "cm")
plotDat = do.call(rbind,performanceCrossTissueSubstType)
plotDat$trainTissue = t2T[plotDat$trainTissue]
plotDat$testTissue = t2T[plotDat$testTissue]
plotDat$substitution = mutClass[plotDat$substitution,1]
ggplot(plotDat, aes(trainTissue,testTissue , fill= AUC)) +
  geom_tile() +
  # scale_fill_viridis_c(limits=colorRange) +
  scale_fill_viridis_c()+
  theme_minimal() +
  facet_wrap(~substitution, nrow = 2) +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1),
        strip.text = element_text(size=12, face="bold"))+
  labs(x = "Train data tissue", y = "Test data tissue") 
ggsave(filename = paste0("fig/CrossTissue/crossTissueSubstitutionType", 
                         plotEnding, ".png"),
       width = 20, height = 15, units = "cm")
ggsave(filename = paste0("fig/CrossTissue/crossTissueSubstitutionType", 
                         plotEnding, ".pdf"),
       width = 20, height = 15, units = "cm")
#####