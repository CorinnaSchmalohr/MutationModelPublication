# preparation #####
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(readxl)
library(ROCR)
library(ranger) 
source("scripts/05_analysis/00_NamesAndColors.R")
source("lib/general_function.R")

tissues = c("brain", "breast", "esophagus", "kidney", "liver", "ovary", "prostate", "skin")
t2T = t2T[tissues]
chroms = sort(paste0("chr", 1:22))
nThread = 32
plotEnding = "_20241114"
#####


# get predictions from WGS model on WEX positions #####
# for each tissue
WGSmodelOnWEX = sapply(tissues, function(tissue){
  print(tissue)
  # load the WEX data
  load(paste0("data/MutTables/exomeTrainData/",
              tissue, "_Muts_mapped_processed.RData")) # dat, datchroms
  # adapt to match WGS predictors
  dat$transcript = 1
  dat$coding = 1
  dat$cancer_expression_1Mb = dat$cancer_expression
  dat$cancer_expression_100kb = dat$cancer_expression
  dat$cancer_expression_10kb = dat$cancer_expression
  dat$cancer_expression_1kb = dat$cancer_expression
  dat$cancer_expression_1bp = dat$cancer_expression
  dat$healthy_expression_1Mb = dat$healthy_expression
  dat$healthy_expression_100kb = dat$healthy_expression
  dat$healthy_expression_10kb = dat$healthy_expression
  dat$healthy_expression_1kb = dat$healthy_expression
  dat$healthy_expression_1bp = dat$healthy_expression

  if(!tissue %in% c("brain", "lung", "skin")){
    dat$normal_expression_1Mb = dat$normal_expression
    dat$normal_expression_100kb = dat$normal_expression
    dat$normal_expression_10kb = dat$normal_expression
    dat$normal_expression_1kb = dat$normal_expression
    dat$normal_expression_1bp = dat$normal_expression
  }
  
  # iterate through chromosomes and get predictions
  predPerChr = sapply(chroms, function(cr){
    cat(cr, " ")
    testData = dat[datchroms != cr,]
    # load the WGS model
    
    load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_", 
                cr, "_forPrediction.RData")) # rf
    # predict on WEX data
    pred = predict(rf,data = testData, num.threads = nThread, verbose = F)
    # store predictions along with real mutation state
    res = data.frame(pred = pred$predictions[,2],  label = testData$mutated)
    return(res)
  }, simplify = F); cat("\n")
  # concatenate chromosomes
  predConcat = do.call(rbind,predPerChr)
  # compute AUROC
  AUC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "auc")
  return(AUC)
})
save(WGSmodelOnWEX, file = "data/Modeling/WGSmodelOnWEX.RData")
gc()
#####


# get predictions from WEX model on WGS positions #####
# for each tissue
WEXmodelOnWGS = sapply(tissues, function(tissue){
  print(tissue)
  # load the WGS data
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
              tissue, "_mapped_processed.RData")) # dat, datchroms
  # adapt to match WEX predictors
  dat$eQTL_pval_1Mb = dat$eQTL_pval_100kb
  dat$eQTL_slope_1Mb = dat$eQTL_slope_100kb
  dat$cancer_expression = dat$cancer_expression_1bp
  dat$healthy_expression = dat$healthy_expression_1bp
  dat$normal_expression = dat$normal_expression_1bp
  dat$effect = 0
  
  # iterate through chromosomes and get predictions
  predPerChr = sapply(chroms, function(cr){
    cat(cr, " ")
    testData = dat[datchroms != cr,]
    # load the WEX model
    load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
                cr, "_forPrediction.RData")) # rf
    # predict on WGS data
    pred = predict(rf,data = testData, num.threads = nThread, verbose = F)
    # store predictions along with realm utation state
    res = data.frame(pred = pred$predictions[,2],  label = testData$mutated)
    return(res)
  }, simplify = F); cat("\n")
  # concatenate chromosomes
  predConcat = do.call(rbind,predPerChr)
  # compute AUROC
  AUC = performance(prediction(predConcat$pred, 
                               predConcat$label), 
                    "auc")
  return(AUC)
}, simplify = F)
save(WEXmodelOnWGS, file = "data/Modeling/WEXmodelOnWGS.RData")
#####


print("done")

# visualize #####
load("data/Modeling/WGSmodelOnWEX.RData")
WGSmodelOnWEX = sapply(WGSmodelOnWEX, function(x){x@y.values[[1]]})
load("data/Modeling/WEXmodelOnWGS.RData")
WEXmodelOnWGS = sapply(WEXmodelOnWGS, function(x){x@y.values[[1]]})
# load WEX self-performance
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
WEXself = sapply(ROC_PR_RF_concat[tissues], function(x){x$auc@y.values[[1]]})
# load WGS self-performance
load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")
WGSself = sapply(ROC_PR_RF_concat[tissues], function(x){x$auc@y.values[[1]]})

# visualize WGS model on WEX, vs WEX model on WEX, vs WGS on WGS model, WEX on WEX. 
AUCs = rbind(WGSmodelOnWEX, WEXself, WEXmodelOnWGS,  WGSself)
pdfAndPng(file= paste0("fig/WGSonWEXandWEXonWGS", plotEnding), 
          width = 10, height = 4, pngArgs = list(pointsize = 20), 
          pdfArgs = list(pointsize=15), expr = expression({
            par(mar = c(3,4,3,1))
            barplot(AUCs-0.5, beside = T, las = 1, names.arg = t2T[colnames(AUCs)],ylab = "AUC",
                    legend.text = c("WGS model on WEX data","WEX model on WEX data",
                                    "WEX model on WGS data","WGS model on WGS data"), 
                    args.legend = list(x = 35, y = 0.22, inset = 0.01,xpd = NA),
                    col = palette.colors(n = 4,"Paired"), yaxt = "n")
            axis(2, at = axTicks(2), labels = c(0,axTicks(2)[-1]+0.5), las = 1)
            plotrix::axis.break(axis=2, breakpos = 0.01, style = "gap")
          }))

#####