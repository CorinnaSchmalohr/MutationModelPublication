.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(readxl)
library(ROCR)
library(ranger) 
source("scripts/05_analysis/00_NamesAndColors.R")
tissues = c("brain", "breast", "esophagus", "kidney", "liver", "ovary", "prostate", "skin")
t2T = t2T[tissues]
chroms = sort(paste0("chr", 1:22))

# get predictor order ######
tab = read_xlsx("data/rawdata/dataMappingAlltissues_WGS.xlsx",
                sheet="allTissues", col_names=T)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000,
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
tab$NA. = NULL
tab[tab == "NA"] = NA
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
predictorOrder = tab[,1:3] # Group, Name, abbreviation
p2P = setNames(predictorOrder$Name, predictorOrder$abbreviation)
predictorGroups = split(predictorOrder$abbreviation, predictorOrder$Group)
rm(tab, predictorOrder)
#####

# # get mutation rate and nPositions, percTP, correlation etc. (dataInfos) #####
# print("dataInfos")
# dataInfos = sapply(tissues, function(tissue){
#   print(tissue)
#   load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
#               tissue, "_mapped_processed.RData"))
#   chrCount = sapply(table(datchroms), as.integer)
#   chrPerc = sapply(split(dat$mutated, datchroms),function(x){mean(x==1)})
#   cors = cor(dat[sapply(dat, is.numeric)], use = "pair")
#   return(list(nMuts = nrow(dat),
#               percTP = mean(dat$mutated == 1),
#               nMutsPerChr = chrCount,
#               percTPperChr = chrPerc,
#               cors = cors))
# }, simplify = F)
# save(dataInfos, file = "data/processedData/WholeGenomedataInfos.RData")
# #####
# 
# # load predictions for each model and each tissue #####
# print("predictions")
# # random forest
# predPerTissueRF = sapply(tissues, function(tissue){
#   print(tissue)
#   load(paste0("data/Modeling/WholeGenomeData/RF/", tissue,
#               "_predictions.RData"))
#   return(predictions)
# }, simplify=F)
# save(predPerTissueRF,
#      file = "data/Modeling/WholeGenomeData/RF/predPerTissueRF.RData")
# # load sig. GLM predictions for each tissue
predPerTissueGLMsig = sapply(tissues, function(tissue){
  load(file = paste0("data/Modeling/WholeGenomeData/GLM/", tissue,
                     "_predictions_sig.RData"))
  return(predictions)
}, simplify=F)
save(predPerTissueGLMsig,
     file = "data/Modeling/WholeGenomeData/GLM/predPerTissueGLM_sig.RData")
# # LASSO
# # predPerTissueLasso = sapply(tissues[1:7], function(tissue){
# #   load(file = paste0("data/Modeling/WholeGenomeData/Lasso/",
# #                      tissue, "_predictions_sig.RData"))
# #   return(predictions)
# # }, simplify=F)
# # save(predPerTissueLasso,
# #      file = "data/Modeling/WholeGenomeData/Lasso/predPerTissueLasso.RData")
# #####
# 
# 
# # compute ROC, PR, and AUROC for each chromosome  #####
# print("ROC PR")
# # RF
# ROC_PR_RF_perChr = sapply(names(predPerTissueRF), function(tissue){ #iterate through tissues
#   pred = predPerTissueRF[[tissue]]
#   # get performances
#   res = lapply(pred, function(x){ #iterate through chromosomes
#     perf = prediction(x$pred, x$label)
#     roc = performance(perf, "tpr", "fpr")
#     pr = performance(perf,"prec", "rec")
#     auc = performance(perf,"auc")@y.values[[1]]
#     return(list(roc = roc, pr = pr, auc = auc))
#   })
#   return(res)
# }, simplify=F)
# save(ROC_PR_RF_perChr,
#      file = "data/Modeling/WholeGenomeData/RF/ROC_PR_RF_perChr.RData")
# # glm
ROC_PR_glm_perChr_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  pred = predPerTissueGLMsig[[tissue]]
  # get performances
  res = lapply(pred, function(x){ #iterate through chromosomes
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  })
  return(res)
}, simplify=F)
save(ROC_PR_glm_perChr_sig,
     file = "data/Modeling/WholeGenomeData/GLM/ROC_PR_glm_perChr_sig.RData")
# # lasso
# # ROC_PR_lasso_perChr = sapply(names(predPerTissueLasso), function(tissue){
# #   pred = predPerTissueLasso[[tissue]]
# #   # get performances
# #   res = lapply(pred, function(x){ #iterate through chromosomes
# #     perf = prediction(x$pred, x$label)
# #     roc = performance(perf, "tpr", "fpr")
# #     pr = performance(perf,"prec", "rec")
# #     auc = performance(perf,"auc")@y.values[[1]]
# #     return(list(roc = roc, pr = pr, auc = auc))
# #   })
# #   return(res)
# # }, simplify=F)
# # save(ROC_PR_lasso_perChr,
# #      file = "data/Modeling/WholeGenomeData/Lasso/ROC_PR_lasso_perChr.RData")
# #####
# 
# 
# # compute ROC, PR, and AUROC for all chromosomes concatenated #####
# print("ROC PR concat")
# ROC_PR_RF_concat = sapply(names(predPerTissueRF), function(tissue){
#   predConcat = do.call(rbind,predPerTissueRF[[tissue]])
#   ROC_rf = performance(prediction(predConcat$pred,
#                                   predConcat$label),
#                        "tpr", "fpr")
#   AUC_rf = performance(prediction(predConcat$pred,
#                                   predConcat$label),
#                        "auc")
#   PR_rf = performance(prediction(predConcat$pred,
#                                  predConcat$label),
#                       "prec", "rec")
#   return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
# }, simplify=F)
# save(ROC_PR_RF_concat,
#      file = "data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")
# #  glm
ROC_PR_glm_concat_sig = sapply(names(predPerTissueGLMsig), function(tissue){
  predConcat = do.call(rbind,predPerTissueGLMsig[[tissue]])
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
save(ROC_PR_glm_concat_sig,
     file = "data/Modeling/WholeGenomeData/GLM/ROC_PR_glm_concat_sig.RData")
# # lasso
# # ROC_PR_lasso_concat = sapply(names(predPerTissueLasso), function(tissue){
# #   predConcat = do.call(rbind,predPerTissueLasso[[tissue]])
# #   ROC_rf = performance(prediction(predConcat$pred,
# #                                   predConcat$label),
# #                        "tpr", "fpr")
# #   AUC_rf = performance(prediction(predConcat$pred,
# #                                   predConcat$label),
# #                        "auc")
# #   PR_rf = performance(prediction(predConcat$pred,
# #                                  predConcat$label),
# #                       "prec", "rec")
# #   return(list(roc = ROC_rf, pr = PR_rf, auc = AUC_rf))
# # }, simplify=F)
# # save(ROC_PR_lasso_concat,
# #      file = "data/Modeling/WholeGenomeData/Lasso/ROC_PR_lasso_concat.RData")
# #####
# 
# 
# # compute concat train performance of RF, GLM and LASSO #####
print("ROCPR train")
print("RF")
predPerTissueRF_train = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
              tissue, "_mapped_processed.RData"))
  chroms = unique(datchroms)
  predConcat = sapply(chroms, function(cr){
    cat(cr,' ')
    load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_",
                cr, "_forPrediction.RData"))
    temp = data.frame(pred = rf$predictions[,2],  label = dat[datchroms != cr,]$mutated)
  }, simplify = F); cat('\n')
  predConcat = do.call(rbind,predConcat)
  return(predConcat)
}, simplify = F)
print("now ROC")
ROC_PR_RF_train = sapply(tissues, function(tissue){
  print(tissue)
  ROC = performance(prediction(predPerTissueRF_train[[tissue]]$pred,
                               predPerTissueRF_train[[tissue]]$label),
                    "tpr", "fpr")
  AUC = performance(prediction(predPerTissueRF_train[[tissue]]$pred,
                               predPerTissueRF_train[[tissue]]$label),
                    "auc")
  PR = performance(prediction(predPerTissueRF_train[[tissue]]$pred,
                              predPerTissueRF_train[[tissue]]$label),
                   "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify = F)
save(ROC_PR_RF_train,
     file = "data/Modeling/WholeGenomeData/RF/ROC_PR_RF_train.RData")
print("GLM")
ROC_PR_glm_sig_train = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
              tissue, "_mapped_processed.RData"))
  chroms = unique(datchroms)
  pred = sapply(chroms, function(cr){
    trainData = dat[datchroms != cr,]
    load(paste0("data/Modeling/WholeGenomeData/GLM/",
                tissue, "_", cr, "_sig.RData"))
    temp = data.frame(pred = logR$fitted.values,  label = trainData$mutated)
  }, simplify = F)
  predConcat = do.call(rbind,pred)
  ROC = performance(prediction(predConcat$pred,
                               predConcat$label),
                    "tpr", "fpr")
  AUC = performance(prediction(predConcat$pred,
                               predConcat$label),
                    "auc")
  PR = performance(prediction(predConcat$pred,
                              predConcat$label),
                   "prec", "rec")
  return(list(roc = ROC, pr = PR, auc = AUC))
}, simplify = F)
save(ROC_PR_glm_sig_train,
     file = "data/Modeling/WholeGenomeData/RF/ROC_PR_glm_sig_train.RData")
# print("LASSO")
# ROC_PR_lasso_sig_train = sapply(tissues[1:7], function(tissue){
#   print(tissue)
#   load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
#               tissue, "_mapped_processed.RData"))
#   chroms = unique(datchroms)
#   pred = sapply(chroms, function(cr){
#     trainData = dat[datchroms != cr,]
#     load(paste0("data/Modeling/WholeGenomeData/Lasso/",
#                 tissue, "_", cr, "_sig.RData"))
#     temp = data.frame(pred = logR$fitted.values,  label = trainData$mutated)
#   }, simplify = F)
#   predConcat = do.call(rbind,pred)
#   ROC = performance(prediction(predConcat$pred,
#                                predConcat$label),
#                     "tpr", "fpr")
#   AUC = performance(prediction(predConcat$pred,
#                                predConcat$label),
#                     "auc")
#   PR = performance(prediction(predConcat$pred,
#                               predConcat$label),
#                    "prec", "rec")
#   return(list(roc = ROC, pr = PR, auc = AUC))
# }, simplify = F)
# save(ROC_PR_lasso_sig_train,
#      file = "data/Modeling/WholeGenomeData/RF/ROC_PR_lasso_sig_train.RData")
#####


# get predictor importances #####
# print("importances")
# # rf
# print("rf")
rf_gini = sapply(tissues, function(tissue){
  load(paste0("data/Modeling/WholeGenomeData/RF/", tissue,
              "_importances_gini.RData"))
  names(imp) = names(chrCols)
  imp = as.data.frame(imp)
  temp = imp[names(p2P),]
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
save(rf_gini, file = "data/Modeling/WholeGenomeData/RF/RF_imps.RData")
# glm significant coefficients
print("glm")
glm_imps_sig = sapply(tissues, function(tissue){
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/WholeGenomeData/GLM/",
                tissue, "_", cr, "_sig.RData"))
    logR$coefficients[names(p2P)]
  })
  rownames(temp) = names(p2P)
  return(temp)
}, simplify=F)
# glm significant pvals
print("glmpvals")
glm_pvals_sig = sapply(tissues, function(tissue){
  print(tissue)
  temp = sapply(names(chrCols),function(cr){
    load(paste0("data/Modeling/WholeGenomeData/GLM/",
                tissue, "_", cr, "_sig.RData"))
    pvals = coef(summary(logR))[-1,4][names(p2P)]
  })
  temp[temp == 0] = 2e-16
  rownames(temp) = names(p2P)

  return(temp)
}, simplify=F)
save(glm_imps_sig, glm_pvals_sig,
     file = "data/Modeling/WholeGenomeData/GLM/GLM_impsAndPvals_sig.RData")
# lasso
# print("lasso")
# lasso_stability = sapply(tissues[1:7], function(tissue){
#   temp = sapply(names(chrCols),function(cr){
#     load(paste0("data/Modeling/WholeGenomeData/Lasso/",
#                 tissue, "_", cr, ".RData")) # sp and stab
#     sp$x[,stab$lpos][names(p2P)]
#   })
#   rownames(temp) = names(p2P)
#   return(temp)
# }, simplify=F)
# print("lasso imp")
# lasso_imp = sapply(tissues[1:7], function(tissue){
#   temp = sapply(names(chrCols),function(cr){
#     load(paste0("data/Modeling/WholeGenomeData/Lasso/",
#                 tissue, "_", cr, "_sig.RData")) # logR, sigFeatures
#     logR$coefficients[names(p2P)]
#   })
#   rownames(temp) = names(p2P)
#   return(temp)
# }, simplify=F)
# save(lasso_stability, lasso_imp,
#      file = "data/Modeling/WholeGenomeData/Lasso/Lasso_impAndStab.RData")
#####

# nPerTissue #####
print("nPerTissue")
nPerTissue = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/WholeGenomeData/WGSMuts_", tissue, "_nPerSample.RData"))
  return(counts)
}, simplify = F)
save(nPerTissue, file = "data/MutTables/WholeGenomeData/nPerTissue.RData")
#####


