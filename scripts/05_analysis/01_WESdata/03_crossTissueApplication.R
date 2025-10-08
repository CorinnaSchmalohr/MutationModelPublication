library(ranger)
library(ROCR)
library(ggplot2)
source("scripts/05_analysis/00_NamesAndColors.R")
cr= "chr1"
load("data/processedData/dataInfos.RData")
dir.create("data/Modeling/exomeTrainData/CrossTissue", showWarnings = F)
dir.create("fig/CrossTissue", showWarnings = F)
plotEnding = "_20230902"
colorRange = c(0.45,0.7) 


# apply models to data from other tissues and get performance #####
RFperformanceCrossTissue = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue2predict, "_Muts_mapped_processed.RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/Modeling/exomeTrainData/RF/", predictingTissue, "_", 
                cr, "_forPrediction.RData"))
    load(paste0("data/MutTables/exomeTrainData/", 
                predictingTissue, "_Muts_mapped_processed.RData"))
    predictors = rf$forest$independent.variable.names
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){ # Add mean value of predictors missing for the tissue to predict on (from predicting tissue)
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(rf,data=testData,type="response", num.threads=10)
    temp = prediction(yHat$predictions[,2], truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
PerfsSelf = sapply(tissues, function(tissue){
  RFperformanceCrossTissue[[tissue]][[tissue]]
}, simplify = F)
save(RFperformanceCrossTissue, PerfsSelf, 
     file = "data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissue.RData")
######

# visualize results ######
# load("data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissue.RData")
# plot ROC
png(paste0("fig/CrossTissue/modelOnOtherTissues_ROC", plotEnding,".png"),
    width=1200, height=1200, pointsize=27)
par(mfrow = c(3,4), mar = c(3,3,1,1))
for(predictingTissue in tissues){
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR", 
       main = paste0(t2T[predictingTissue],
                     " model"), las = 1, mgp = c(3,1,0))
  abline(0,1)
  for(tissue2predict in tissues){
    temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
    plot(temp$roc, add = T, col = tissueCols[tissue2predict], lwd = 2)
    if(tissue2predict != predictingTissue){
      temp2 = PerfsSelf[[tissue2predict]]
      plot(temp2$roc, add = T, col = tissueCols[tissue2predict], 
           lty = 2, lwd = 2)
    }
  }
}
plot.new()
legend("center",legend=t2T[names(tissueCols)], col = tissueCols,cex = 1.2,
       lty = 1, bty = "n", lwd = 2)
plot.new()
legend("center", bty = "n",cex = 1.2,
       legend = c(paste0('"foreign"', " model"), 
                  "tissue-specific model"), lty = 1:2, lwd = 2)
dev.off()


# plot PR
png(paste0("fig/CrossTissue/modelOnOtherTissues_PR", plotEnding,".png"),
    width=1200, height=1200, pointsize=27)
par(mfrow = c(3,4), mar = c(4,4,1,1))
for(predictingTissue in tissues){
  plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall", ylab = "Precision", 
       main = paste0(t2T[predictingTissue],
                     " model"), las = 1, mgp = c(2.5,1,0))
  for(tissue2predict in tissues){
    temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]
    plot(temp$pr, add = T, col = tissueCols[tissue2predict], lwd = 2)
    if(tissue2predict != predictingTissue){
      temp2 = PerfsSelf[[tissue2predict]]
      plot(temp2$pr, add = T, col = tissueCols[tissue2predict],
           lty = 2, lwd = 2)
    }
  }
}
plot.new()
legend("center",legend=t2T[names(tissueCols)], col = tissueCols,cex = 1.2,
       lty = 1, bty = "n", lwd = 2)
plot.new()
legend("center", bty = "n", cex = 1.2,
       legend = c(paste0('"foreign"', " model"),
                  "tissue-specific model"), lty = 1:2, lwd = 2)
dev.off()


# image of ROCs
ROCs = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = RFperformanceCrossTissue[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs$AUC = as.numeric(ROCs$AUC)
ROCs$rank <- NA
for(tissue in tissues){
  ROCs[ROCs$testedOn == tissue,]$rank[order(-ROCs[ROCs$testedOn == tissue,]$AUC)] = 1:length(tissues)
}
ggplot(ROCs, aes(trainedOn, testedOn, fill= AUC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_AUCheatmap", plotEnding,".png"), 
       width = 12, height = 12, units = "cm")
ggplot(ROCs, aes(trainedOn, testedOn, fill= AUC)) +
  geom_tile() +
  coord_equal() + # for square shapes
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "points"))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels= t2T) + 
  scale_y_discrete(labels = t2T)
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_AUCheatmap_simple", plotEnding,".png")) 
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_AUCheatmap_simple", plotEnding,".pdf"))
#####


# apply models to subsampled data from other tissues and get performance #####
RFperformanceCrossTissue_subsampled = sapply(tissues, function(tissue2predict){
  cat("data: ", tissue2predict, "\n")
  load(paste0("data/MutTables/exomeTrainData_subsampled/", 
              tissue2predict, "_Muts_mapped_processed.RData"))
  tempData = dat[datchroms == cr,]
  truePreds = tempData$mutated
  cat("model: ")
  perfs = sapply(tissues, function(predictingTissue){
    cat(predictingTissue, " ")
    load(paste0("data/Modeling/exomeTrainData_subsampled/RF/", predictingTissue, "_", 
                cr, "_forPrediction.RData"))
    load(paste0("data/MutTables/exomeTrainData_subsampled/", 
                predictingTissue, "_Muts_mapped_processed.RData"))
    predictors = rf$forest$independent.variable.names
    testData = tempData[,colnames(tempData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(rf,data=testData,type="response", num.threads=10)
    temp = prediction(yHat$predictions[,2], truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(RFperformanceCrossTissue_subsampled,
     file = "data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissue_subsampled.RData")
#####


#  heatmap of downsampled results ######
# load("data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissue_subsampled.RData")
ROCs_subsampled = data.frame(do.call(rbind,lapply(tissues, function(tissue2predict){
  t(sapply(tissues, function(predictingTissue){
    temp = RFperformanceCrossTissue_subsampled[[tissue2predict]][[predictingTissue]]$auc@y.values[[1]]
    return(c("trainedOn" = predictingTissue, 
             "testedOn" = tissue2predict, "AUC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
ROCs_subsampled$AUC = as.numeric(ROCs_subsampled$AUC)
ROCs_subsampled$rank <- NA
for(tissue in tissues){
  ROCs_subsampled[ROCs_subsampled$testedOn == tissue,]$rank[order(-ROCs_subsampled[ROCs_subsampled$testedOn == tissue,]$AUC)] = 1:length(tissues)
}
ggplot(ROCs_subsampled, aes(trainedOn, testedOn, fill= AUC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_subsampled_AUCheatmap", plotEnding,".png"), 
       width = 12, height = 12, units = "cm")
ggplot(ROCs_subsampled, aes(trainedOn, testedOn, fill= AUC)) +
  geom_tile() +
  coord_equal() + # for square shapes
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "points"))+
  labs(y = "Tested on ", x = "Trained on") +  
  scale_x_discrete(labels= t2T) + 
  scale_y_discrete(labels = t2T)
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_subsampled_AUCheatmap_simple", plotEnding,".png")) 
ggsave(filename=paste0("fig/CrossTissue/modelOnOtherTissues_subsampled_AUCheatmap_simple", plotEnding,".pdf"))

#####

# exchange predictors
# variant 1: e.g. test lung model on lung mutations with kidney predictors #####
# muttissue: model, predictors, predtissue: predictors
CrossTissuePreds = sapply(tissues, function(muttissue){
  cat("muts and model from ", muttissue, "\n")
  # load RF
  load(paste0("data/Modeling/exomeTrainData/RF/", muttissue, "_", 
              cr, "_forPrediction.RData"))
  predictors = rf$forest$independent.variable.names
  
  cat("predictors from ")
  perfs = sapply(tissues, function(predtissue){
    cat(predtissue, " ")
    if(predtissue == muttissue){
      load(paste0("data/MutTables/exomeTrainData/", 
                  muttissue, "_Muts_mapped_processed.RData"))
    }else{
      load(paste0("data/MutTables/exomeTrainDataCrossTissue/", 
                  muttissue, "With",predtissue,"Preds_processed", ".RData"))
    }
    
    testData = dat[datchroms == cr,]
    truePreds = testData$mutated
    # add dummy variables for predictors present only in predtissue
    testData = testData[,colnames(testData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    load(paste0("data/MutTables/exomeTrainData/", 
                muttissue, "_Muts_mapped_processed.RData"))
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(rf,data=testData,type="response", num.threads=10)
    temp = prediction(yHat$predictions[,2], truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(CrossTissuePreds, 
     file = "data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissuePredsExchanged1.RData")
#####


# variant 2: e.g. test lung model on kidney mutations with lung preds #####
print("variant 2")
CrossTissueMuts = sapply(tissues, function(muttissue){
  cat("preds and model from ", muttissue, "\n")
  # load RF
  load(paste0("data/Modeling/exomeTrainData/RF/", muttissue, "_", 
              cr, "_forPrediction.RData"))
  predictors = rf$forest$independent.variable.names
  
  cat("muts from ")
  perfs = sapply(tissues, function(predtissue){
    cat(predtissue, " ")
    if(predtissue == muttissue){
      load(paste0("data/MutTables/exomeTrainData/", 
                  predtissue, "_Muts_mapped_processed.RData"))
    }else{
      load(paste0("data/MutTables/exomeTrainDataCrossTissue/", 
                  predtissue, "With",muttissue,"Preds_processed", ".RData"))
    }
    testData = dat[datchroms == cr,]
    truePreds = testData$mutated
    # add dummy variables for predictors present only in predtissue
    testData = testData[,colnames(testData) %in% predictors]
    toAdd = predictors[!predictors %in% colnames(testData)]
    load(paste0("data/MutTables/exomeTrainData/", 
                muttissue, "_Muts_mapped_processed.RData"))
    for(x in toAdd){
      testData[x] = rep(mean(dat[,x]), nrow(testData))
    }
    yHat = predict(rf,data=testData,type="response", num.threads=10)
    temp = prediction(yHat$predictions[,2], truePreds)
    roc = performance(temp, "tpr", "fpr")
    pr = performance(temp,"prec", "rec")
    auc = performance(temp,"auc")
    perf = list(roc = roc, pr = pr, auc = auc)
    return(perf)
  }, simplify=F)
  cat("\n")
  return(perfs)
}, simplify = F)
save(CrossTissueMuts, 
     file = "data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissuePredsExchanged2.RData")
#####


# visualize results #####
# variant 1
load("data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissuePredsExchanged1.RData")
RocsPreds = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissuePreds[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("TrainedOn" = muttissue, 
             "TestedOn" = predtissue, "AUC" = temp))
    # return(c("predictorSource" = predtissue, 
    #          "MutationAndModelSource" = muttissue, "AUC" = temp))
  }, USE.NAMES= F))  
})), stringsAsFactors=F)
RocsPreds$AUC = as.numeric(RocsPreds$AUC)
RocsPreds$rankColumn <- NA
for(tissue in tissues){
  RocsPreds$rankColumn[RocsPreds$TrainedOn == tissue] = rank(-RocsPreds$AUC[RocsPreds$TrainedOn == tissue])  
}
RocsPreds$rankRow <- NA
for(tissue in tissues){
  RocsPreds$rankRow[RocsPreds$TestedOn == tissue] = rank(-RocsPreds$AUC[RocsPreds$TestedOn == tissue])  
}
ggplot(RocsPreds, aes(y = TestedOn, x = TrainedOn, fill= AUC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
  labs(y = "Test data predictors from this tissue", x = "Model trained on this tissue\n and test data positions from this tissue") +
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rankColumn), color = "white", size = 2.5)
  # geom_text(aes(label = rankRow), color = "white", size = 2.5)
ggsave(filename=paste0("fig/CrossTissue/CrossTissuePredsExchanged", plotEnding,".png"), 
       width = 12, height = 12, units = "cm")
ggplot(RocsPreds, aes(y = TestedOn, x = TrainedOn, fill= AUC)) +
  geom_tile() +
  coord_equal() + # for square shapes
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "points"))+
  labs(y = "Test data predictors", x = "Model and test data positions") +
  scale_x_discrete(labels= t2T) + 
  scale_y_discrete(labels = t2T)
ggsave(filename=paste0("fig/CrossTissue/CrossTissuePredsExchanged_simple", plotEnding,".png")) 
ggsave(filename=paste0("fig/CrossTissue/CrossTissuePredsExchanged_simple", plotEnding,".pdf"))



# variant 2
load("data/Modeling/exomeTrainData/CrossTissue/RFperformanceCrossTissuePredsExchanged2.RData")
RocsMuts = data.frame(do.call(rbind,lapply(tissues, function(muttissue){
  t(sapply(tissues, function(predtissue){
    temp = CrossTissueMuts[[muttissue]][[predtissue]]$auc@y.values[[1]]
    return(c("trainedOn" = muttissue, 
             "testedOn" = predtissue, "AUC" = temp))
  }, USE.NAMES= F))
})), stringsAsFactors=F)
RocsMuts$AUC = as.numeric(RocsMuts$AUC)
RocsMuts$rank <- NA
for(tissue in tissues){
  RocsMuts$rank[RocsMuts$testedOn == tissue] = rank(-RocsMuts$AUC[RocsMuts$testedOn == tissue])  
}
ggplot(RocsMuts, aes(y = testedOn, x = trainedOn, fill= AUC)) +
  geom_tile() +
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1))+
  labs(y = "Test data mutations from this tissue", x = "Model trained on this tissue\n and test data predictors from this tissue") + 
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  geom_text(aes(label = rank), color = "white", size = 2.5)
ggsave(filename=paste0("fig/CrossTissue/CrossTissueMutsExchanged", plotEnding,".png"), 
       width = 12, height = 12, units = "cm")
ggplot(RocsMuts, aes(y = testedOn, x = trainedOn, fill= AUC)) +
  geom_tile() +
  coord_equal() + # for square shapes
  scale_fill_viridis_c(limits=colorRange) +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 46, vjust = 1, hjust = 1),
        axis.text = element_text(size=20),
        axis.title=element_text(size=20),
        legend.text = element_text(size=15),
        legend.title = element_text(size=20),
        legend.key.size = unit(25, "points"))+
  labs(y = "Test data mutations", x = "Model and test data predictors") + 
  scale_x_discrete(labels= t2T) + 
  scale_y_discrete(labels = t2T)
ggsave(filename=paste0("fig/CrossTissue/CrossTissueMutsExchanged_simple", plotEnding,".png")) 
ggsave(filename=paste0("fig/CrossTissue/CrossTissueMutsExchanged_simple", plotEnding,".pdf"))
#####

