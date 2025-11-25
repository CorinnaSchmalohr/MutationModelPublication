source("scripts/05_analysis/00_NamesAndColors.R")
source("lib/general_function.R")

# demonstrate data size differences #####
load("data/MutTables/exomeTrainData_subsampled/nMuts.RData") # nMuts,sampleN
nMuts = nMuts[tissues]
load("data/MutTables/exomeTrainData_subsampledNsamples/nSamps.RData") # nSamps,nSamp

pdfAndPng(file = "fig/exomeData_downsampling_nMutcomparison", 
          width = 12, height = 4, pngArgs = list(pointsize = 20), 
          pdfArgs = list(pointsize=15), 
          expr = expression({
            par(mfrow = c(1,2),mar = c(3,1,1,1), oma = c(0,6,0,0))
            barplot(rev(nMuts), xlab = "n Positions", las = 1, names.arg = rev(t2T[tissues]), 
                    col = rev(tissueCols), mgp =  c(2,0.7,0), horiz = T, xaxt = "n")
            axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = T), mgp =  c(2,0.7,0))
            abline(v = sampleN, lwd = 1.5, lty = 2)
            abline(v=0)
            # par(mar = c(3,0,1,7))
            barplot(rev(nSamps), xlab = "n Samples", las = 1, names.arg = "", 
                    col = rev(tissueCols), mgp =  c(2,0.7,0), horiz = T, xaxt = "n")
            axis(1, at = axTicks(1), labels = format(axTicks(1), scientific = T), mgp =  c(2,0.7,0))
            abline(v = nSamp, lwd = 1.5, lty = 2)
            abline(v=0)
          }))

# dev.off()
# png("fig/exomeData_downsampling_nSampcomparison.png", width = 600, height = 600, pointsize = 14.5)
# par(mar = c(3,6,1,2))
# 
# dev.off()
#####


# load performance for full dataset
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
######


# compute performance for downsampled dataset #####
predPerTissueRF_downsampled = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData_subsampled/RF/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
ROC_PR_RF_concat_downsampled = sapply(names(predPerTissueRF_downsampled), function(tissue){
  predConcat = do.call(rbind,predPerTissueRF_downsampled[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_RF_concat_downsampled, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat_downsampled.RData")
#####


# compute performance for downsampledNsamples dataset #####
predPerTissueRF_downsampledNsamples = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData_subsampledNsamples/RF/", tissue,
              "_predictions.RData"))
  return(predictions)
}, simplify=F)
ROC_PR_RF_concat_downsampledNsamples = sapply(names(predPerTissueRF_downsampledNsamples), function(tissue){
  predConcat = do.call(rbind,predPerTissueRF_downsampledNsamples[[tissue]])
  ROC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "tpr", "fpr")
  AUC_rf = performance(prediction(predConcat$pred, 
                                  predConcat$label), 
                       "auc")
  PR_rf = performance(prediction(predConcat$pred, 
                                 predConcat$label), 
                      "prec", "rec")
  return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
}, simplify=F)
save(ROC_PR_RF_concat_downsampledNsamples, 
     file = "data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat_downsampledNsamples.RData")
#####



# compare performances (AUC) #####
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat_downsampled.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat_downsampledNsamples.RData")
AUCs = sapply(tissues, function(tissue){
  AUCtissue = c(full = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
                downsampled = ROC_PR_RF_concat_downsampled[[tissue]]$auc@y.values[[1]],
                downsampledNsamples = ROC_PR_RF_concat_downsampledNsamples[[tissue]]$auc@y.values[[1]])
  return(AUCtissue)
})
tempCols = rbind(tissueCols, addAlpha(tissueCols, alpha = 0.6), addAlpha(tissueCols, alpha = 0.2))
pdfAndPng(file = "fig/exomeData_downsampling_AUCcomparison", 
          width = 10, height = 4, pngArgs = list(pointsize = 14), pdfArgs = list(pointsize=10),
          expr = expression({
            par(mar = c(4,4,1,1))
            barplot(AUCs-0.5, beside = T, las = 1,  ylab = "AUC", yaxt = "n",
                    names.arg = t2T[tissues], 
                    legend.text = c("Full data", "Downsampled mutations", "Downsampled samples"),
                    args.legend = list(x= "topleft", inset = 0.05, fill = c("black", "grey55", "grey85")),
                    col = tempCols)
            axis(2, at = axTicks(2), labels = c(0,axTicks(2)[-1]+0.5), las = 1)
            plotrix::axis.break(axis = 2, breakpos = 0.01, style = "gap")
          }))
# 
# dev.off()
#####