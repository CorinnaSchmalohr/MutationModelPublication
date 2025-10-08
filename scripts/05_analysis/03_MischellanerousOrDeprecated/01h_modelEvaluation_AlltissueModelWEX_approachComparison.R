.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ggplot2)
library(ROCR)
library(ranger) 
source("lib/general_function.R")
source("scripts/05_analysis/00_NamesAndColors.R")
plotEnding = "_20241024"
tissueCols = c(tissueCols, 
               "TissueCombination" = "purple4",
               "TissueCombination2" = "firebrick4", 
               "allTissues" = "darkgreen")
t2T = c(t2T, "TissueCombination" = "Combined data",
        "TissueCombination2" = "Combined reduced predictors", 
        "allTissues" = "Alltissue data")
approaches = c("TissueCombination", "TissueCombination2", "allTissues")

nThread = 14
# 
# # Mutation overview #####
# print("mutation overview")
# # prepare counts of mutation types
# bases = c("A", "C", "G", "T")
# trimClass = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
#                                 c(rep("C", 48), rep("T", 48)), ">",
#                                 c(rep(c("A", "G", "T"), each =16), 
#                                   rep(c("A", "C", "G"), each =16)), "]",
#                                 rep(bases, 24))),
#                   paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
#                                 c(rep("G", 48), rep("A", 48)), ">",
#                                 c(rep(c("T", "C", "A"), each =16), 
#                                   rep(c("T", "G", "C"), each =16)), "]",
#                                 rep(rev(bases), 24))))
# rownames(trimClass) = paste(trimClass[,1], trimClass[,2], sep = " and ")
# trimTypes = sapply(c(tissues, "allTissues  "), function(tissue){
#   load(paste0("data/MutTables/exomeTrainData/", 
#               tissue, "_Muts_mapped.RData"))
#   muts = data$muts[data$muts$mutated == 1,]
#   temp = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
#                 substr(muts$context,4,4))
#   apply(trimClass,1, function(x){
#     sum(temp %in% x)
#   })
# })
# 
# # plot
# png(paste0("fig/MutationOverview_AllTissueModel", plotEnding, ".png"),
#     width = 2900, height = 1500, res=200, pointsize=18)
# layout(rbind(1:11))
# # par(mar = c(0.1,0.15,0,0), oma = c(3.3,5,2,0))
# par(mar = c(0,1,2,0), oma = c(3.3,5,0,0))
# # lower part: frequencies of mutation types
# temp =apply(trimTypes, 2,function(x){x/sum(x)})
# dumpVar = sapply(1:ncol(temp), function(i){
#   barplot(rev(temp[,i]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
#           names.arg="",cex.names=0.5, cex.axis = 1.1,main  = c(t2T, "All tissues")[i],
#           col = rep(rev(basesCol), each=16), xlim = c(0,0.16), xaxp = c(0, 0.1, 1))
#   if(i==1){
#     classes = trimClass[rownames(temp),1]
#     names = rev(paste0(substr(classes,1,1), substr(classes,3,3), substr(classes,7,7)))
#     names2 = rev(paste0(" ", substr(classes,3,3), " "))
#     text(x = 0, y = 0.99:96-0.5,  labels=names,family = "mono",
#          xpd = NA, adj=1.1,  cex = 0.5)
#     text(x = 0, y = 0.99:96-0.5,  labels=names2,family = "mono",
#          xpd = NA, adj=1.1, col = rep(rev(basesCol), each=16), font = 2, cex = 0.5)
#     text(x = -0.055, y = 1:6*16-8,
#          labels = c("T>G", "T>C", "T>A" ,"C>T", "C>G", "C>A"),
#          xpd = NA, adj=1, col = rev(basesCol), font = 2)
#     segments(y0=1:6*16-16, y1 = 1:6*16, x0=-0.045, lwd = 4, 
#              col = rev(basesCol), xpd = NA, lend=1, cex = 1.1)
#     title(ylab = "Mutation type", mgp = c(5,1,0), xpd = NA, cex.lab = 1.1)
#   }
#   abline(h=seq(16,80,length.out = 5), lty = 2, lwd = 1.5)
#   abline(h=0, lty = 1, lwd = 2)
#   #abline(h=96, lty = 1, lwd = 1.8)
# })
# axis(3,at =1:length(tissueCols)*2-0.5, labels= t2T[names(tissueCols)], tick=F,
#      mgp = c(2,0.5,0), cex.axis  = 1.1)
# title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
# dev.off()
# rm(trimClass, bases, trimTypes, temp, dumpVar)
# #####
# 
# # compute performances #####
# ROC_PR_RF_perChr= sapply(approaches, function(approach){
#   load(paste0("data/Modeling/exomeTrainData/RF/", approach,
#               "_predictions.RData")) # predictions
#   perChr  = lapply(predictions, function(x){ #iterate through chromosomes
#     perf = prediction(x$pred, x$label)
#     roc = performance(perf, "tpr", "fpr")
#     pr = performance(perf,"prec", "rec")
#     auc = performance(perf,"auc")@y.values[[1]]
#     return(list(roc = roc, pr = pr, auc = auc))
#   })
#   return(perChr)
# }, simplify = F)
# ROC_PR_RF_concatApproaches = sapply(approaches, function(approach){
#   load(paste0("data/Modeling/exomeTrainData/RF/", approach,
#               "_predictions.RData")) # predictions
#   predConcat = do.call(rbind,predictions)
#   perf = prediction(predConcat$pred, predConcat$label)
#   roc = performance(perf, "tpr", "fpr")
#   pr = performance(perf,"prec", "rec")
#   auc = performance(perf,"auc")@y.values[[1]]
#   return(list(roc = roc, pr = pr, auc = auc))
# }, simplify = F)
# #####
# 
# # ROC and PR of CWCV ####
# print("ROC and PR")
# 
# # plot
# png(paste0("fig/modelEvaluation/ROC_PR_CWCV_allTissueModel", plotEnding, ".png"),
#     width = 2900, height = 3100, res=200, pointsize=24)
# layout(rbind(c(1,3,5), c(2,4,6)))
# dumpVar= lapply(approaches, function(approach){
#   # ROC
#   plot(NA, xlim = c(0,1), ylim = c(0,1),
#        mgp = c(2,0.7,0), las = 1,
#        xlab ="FPR", ylab = "TPR", main = approach)
#   abline(0,1, col = "grey", lty = 2)
#   dump = sapply(names(chrCols), function(cr){
#     lines(ROC_PR_RF_perChr[[approach]][[cr]]$roc@x.values[[1]], 
#           ROC_PR_RF_perChr[[approach]][[cr]]$roc@y.values[[1]], 
#           col = rgb(0,0,0,0.3))
#   })
#   lines(ROC_PR_RF_concatApproaches[[approach]]$roc@x.values[[1]], 
#         ROC_PR_RF_concatApproaches[[approach]]$roc@y.values[[1]], 
#         col = tissueCols[approach], lwd = 3)
#   
#   legend("bottomright", lty = 1, lwd = c(1,3), bty = "n",
#          col = c(rgb(0,0,0,0.5), tissueCols[approach]),
#          legend=c("CWCV", "Total"))
#   #  PR
#   plot(NA, xlim = c(0,1), ylim = c(0,1),
#        mgp = c(2,0.8,0), las = 1,xlab = "Recall", ylab = "Precision")
#   dump = sapply(names(chrCols), function(cr){
#     lines(ROC_PR_RF_perChr[[approach]][[cr]]$pr@x.values[[1]], 
#           ROC_PR_RF_perChr[[approach]][[cr]]$pr@y.values[[1]], 
#           col = rgb(0,0,0,0.3))
#   })
#   lines(ROC_PR_RF_concatApproaches[[approach]]$pr@x.values[[1]], 
#         ROC_PR_RF_concatApproaches[[approach]]$pr@y.values[[1]], 
#         col = tissueCols[approach], lwd = 3)
# })
# dev.off()
# #####
# 
# 
# # comparison of tissue performance #####
# print("comparison of tissues")
# load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
# 
# png(paste0("fig/modelEvaluation/compareRFperformance_allTissueModel", plotEnding, ".png"), 
#     width=2400, height=800, pointsize=35)
# par(mfrow = c(1,3), mar = c(4,7.5,0.1,0.1), mgp = c(2.5,1,0))
# # roc
# plot(NA, xlim = c(0,1), ylim = c(0,1), 
#      xlab = "FPR", 
#      ylab = "TPR", las = 1)
# abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
#         ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
#         col = tissueCols[tissue], lwd = 2)
# })
# plotDump = sapply(approaches, function(approach){
#   lines(ROC_PR_RF_concatApproaches[[approach]]$roc@x.values[[1]], 
#         ROC_PR_RF_concatApproaches[[approach]]$roc@y.values[[1]], 
#         col = tissueCols[approach], lwd = 2)
# })
# # pr
# plot(NA, xlim = c(0,1), ylim = c(0,1), 
#      xlab = "Recall", 
#      ylab = "Precision", las = 1)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]], 
#         ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], 
#         col = tissueCols[tissue], lwd = 2)
# })
# plotDump = sapply(approaches, function(approach){
#   lines(ROC_PR_RF_concatApproaches[[approach]]$pr@x.values[[1]], 
#         ROC_PR_RF_concatApproaches[[approach]]$pr@y.values[[1]], 
#         col = tissueCols[approach], lwd = 2)
# })
# legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)], lwd = 2)
# # AUCs
# AUCs = c(sapply(tissues, function(tissue){
#   ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
# }), sapply(approaches, function(approach){
#   ROC_PR_RF_concatApproaches[[approach]]$auc
# }))
# temp = barplot(AUCs, col = tissueCols, las = 1, ylab = "AUC", names.arg = "")
# text(x = temp, y = -0.02, labels = t2T,srt = 45, c(1,0), 
#      xpd = T)
# dev.off()
# 
# #####
# 
# 
# 
# # comparison of predictor importance #####
# print("comparison of predictor importances")
# rf_imps = do.call(rbind,sapply(names(tissueCols), function(tissue){
#   load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
#               "finalModel.RData")) # rf, importance
#   res = data.frame("tissue" = tissue, 
#                    "predictor" = names(p2P), 
#                    "gini" = importance[names(p2P)],
#                    "gini_scaled" = scale(importance[names(p2P)], center = F))
#   return(res)
# }, simplify=F))
# rf_imps$predictor = factor(rf_imps$predictor, levels = names(p2P))
# ggplot(rf_imps, aes(x = tissue, y = predictor, fill = gini_scaled)) + 
#   geom_raster() +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") +
#   labs(y = "Predictor", x = "Tissue") + 
#   labs(fill = "Scaled Gini\nImportance")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2P)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
#         axis.title=element_text(size=14,face="bold")) 
# ggsave(paste0("fig/modelEvaluation/RFginiScaled_allTissueModel", plotEnding, ".png"), 
#        height=8, width=6)
# rm(rf_imps)
# #####


# apply to other tissues #####
print("apply to other tissues")
crossTissuePerformance = sapply(approaches, function(approach){
  print(approach)
  if(approach == "allTissues"){
    load(paste0("data/MutTables/exomeTrainData/", "allTissues", "_Muts_mapped_processed.RData"))
    trainMeans = sapply(dat[,colnames(dat) != "mutated"], mean, na.rm = T); rm(dat, datchroms)
    predictors = names(trainMeans)
  }
  perf = sapply(tissues, function(tissue){
    print(tissue)
    load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped_processed.RData"))
    # match predictors
    if(approach == "allTissues"){
      toAdd = predictors[!predictors %in% colnames(dat)]
      for(x in toAdd){
        dat[x] = rep(trainMeans[x], nrow(dat))
      }
      dat = dat[,c(names(trainMeans), "mutated")]}
    # iterate through chromosomes
    predPerChr = sapply(names(chrCols), function(cr){
      # load general model
      if(approach == "allTissues"){
        load(paste0("data/Modeling/exomeTrainData/RF/", "alltissues", "_", cr, "_forPrediction.RData"))
      } else{
        load(paste0("data/Modeling/exomeTrainData/RF/", approach, "_", cr, "_forPrediction.RData"))
      }
      # predict on data using model
      testDat = dat[datchroms == cr,]
      yHat = predict(rf,data=testDat,type="response", num.threads=nThread)
      # store predictions + true value
      res = data.frame(prediction = yHat$predictions[,2], 
                       trueVal = as.integer(testDat$mutated)-1)
      return(res)
    }, simplify = F)
    # concatenate results
    preds = do.call(rbind, predPerChr)
    # compute performance
    temp = prediction(predictions = preds$prediction, labels = preds$trueVal)
    auc = performance(temp,"auc")
    return(auc@y.values[[1]])
  })
}, simplify = F)

save(crossTissuePerformance, 
     file = "data/Modeling/exomeTrainData/CrossTissue/allTissue_crossTissuePerformance.RData")
######


#visualize cross-tissue performance #####
load("data/Modeling/exomeTrainData/CrossTissue/allTissue_crossTissuePerformance.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
tissueSpecificPerformance = sapply(tissues, function(tissue){
  ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
})
AUCs = rbind(tissueSpecificPerformance, do.call(rbind,crossTissuePerformance))
png(paste0("fig/modelEvaluation/allTissueModelCrossTissuePerformance", plotEnding, ".png"),
    width=800, height=400, pointsize=12)
par(oma = c(0,0,5,0), mar = c(3,4,0,0))
temp = barplot(AUCs, beside = T,  las = 1, ylab = "AUC",
        names.arg = t2T[tissues], col = rep(tissueCols[1:10], each = 4),
        density = c(-1, 10,20,10), angle = c(-1,0,85,-45), legend.text = F, bg = "grey")
legend(x = 0,y = max(AUCs)+0.16, xjust = 0, yjust = 1,
       legend = c("tissue-specific model","combination of tissue-specific data",
                  "combination with reduced predictors", "allTissue model (non-specific predictors)"),
       density = c(-1, 30,30,20), angle = c(-1,0,85,-45),fill = "black", xpd = NA, inset =0)
dev.off()
#####


