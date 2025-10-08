#library(ranger)
library(ROCR)
library(ggplot2)
library(readxl)
library(gridExtra)
library(RColorBrewer)

source("scripts/05_analysis/00_NamesAndColors.R")
set.seed(123)
dir.create("fig/featureExchange/", showWarnings = F)

tissue_model = tissues[7]
tissue_feature = tissues[5]

# list of significant tissuespecific features of the model 
tab = read_xlsx("data/rawdata/dataMapping.xlsx", 
                sheet=tissue_model, col_names=T)
tissuespec_features = tab$abbreviation[tab$sourceTissue == tissue_model]
#f = features[2]
load(paste0("data/rdata/GLMmodel/", tissue_model,"_pvals.RData"))
sig_features = rowMeans(pvals)
sig_features = names(sig_features[sig_features <= 0.05])


# Load data from which the features to be replaced should come from 
load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",
       tissue_model, "With",tissue_feature,"Preds_processed", ".RData"))
f_dat = dat   # f_dat = data from where to take the features
rm(dat)

# Only check for features with given tissuespecific data 
features = sig_features[sig_features %in% tissuespec_features]
features = features[features %in% colnames(f_dat)]
#missing_f = setdiff(sig_features, features)

# Load tissue data on which the model was trained on  --> dat
load(paste0("data/procData/traindata/traindata_processed_",
            tissue_model, ".RData"))  
chroms = unique(datchroms)

# Iterate through all tissue spec. features #
sapply(features, function(f){
  print(f)
  repl_feature = f_dat[f]
  old_feature = dat[f]
  
  ### Data with only sig. features 
  p_dat = dat   # p_dat = data to make the predictions on
  p_dat = cbind(p_dat[sig_features], mutated = p_dat$mutated)
  # Permute feature of interest 
  p_dat[f] = dat[sample(1:nrow(dat), nrow(dat)), f]

  # Get predictions on data with permuted feature
  print("making predictions on permuted feature...")
  predictions_permuted = lapply(chroms, function(cr){
    load(paste0("data/rdata/GLMmodel/", tissue_model, "_", cr, "_sig.RData")) 
    testData = p_dat[datchroms == cr,]
    yhat = predict(logR, newdata = testData, type = "response")
    temp = data.frame(pred  = yhat,label = testData$mutated)
    return(temp)
  })
  names(predictions_permuted) = chroms
  
  
  
  # Replace feature of interest
  r_dat = cbind(dat[sig_features], mutated = dat$mutated)
  r_dat["diff"] = abs(old_feature - repl_feature) # add difference between features as column
  r_dat[f] = repl_feature

  ### Get predictions on data with exchanged feature ###
  # predict on held-out chromosomes
  print("making predictions...")
  predictions_replF = lapply(chroms, function(cr){
    load(paste0("data/rdata/GLMmodel/", tissue_model, "_", cr, "_sig.RData")) 
    testData = r_dat[datchroms == cr, names(r_dat) != "diff"] # Ignore differnce for prediction 
    yhat = predict(logR, newdata = testData, type = "response")
    temp = data.frame(pred = yhat,  label = testData$mutated, 
                      diff = r_dat[datchroms == cr, names(r_dat) == "diff"]) # Add difference to the output
    return(temp)
  })
  names(predictions_replF) = chroms

  ### General comparison ###
  # compute ROC, PR, and AUROC for all chromosomes concatenated 
  predConcat_replF = do.call(rbind,predictions_replF)
  predConcat_permuted = do.call(rbind,predictions_permuted)
    
  AUC_replF = performance(prediction(predConcat_replF$pred, predConcat_replF$label), "auc")
  AUC_replF = AUC_replF@y.values[[1]]
  
  AUC_permuted = performance(prediction(predConcat_permuted$pred, predConcat_permuted$label), "auc")
  AUC_permuted = AUC_permuted@y.values[[1]]
  
  
  ## Add concatenated predictions of non-changed data ##
  load("data/rdata/GLMmodel/predPerTissueGLM_sig.RData")
  predConcat = do.call(rbind,predPerTissueGLMsig[[tissue_model]])
  
  # Compare to non-changed data ROC and PR
  load(file = "data/rdata/GLMmodel/ROC_PR_glm_concat_sig.RData")
  
  # Combine data for plotting
  combPredsConcat = cbind(predConcat, pred_replF = predConcat_replF$pred, Fdiff = predConcat_replF$diff)
  combPredsConcat["PredDiff"] = combPredsConcat$pred - combPredsConcat$pred_replF
  # Impact on positions with the most difference 
  combPredsConcat = combPredsConcat[order(combPredsConcat$Fdiff, decreasing = T),]
  
  # Sep data in n parts and compare the results of each part
  n = 6
  seg = factor(rep(x = 1:n, each = ceiling(nrow(dat)/6)))
  combPredsConcat["seg"] = seg[1:nrow(combPredsConcat)]

  res_sep = sapply(1:n, function(i){
    data = combPredsConcat[combPredsConcat["seg"] == i,]
    AUC_replF = performance(prediction(data$pred_replF, data$label), "auc")
    AUC = performance(prediction(data$pred, data$label), "auc")
    AUC_diff = AUC_replF@y.values[[1]] - AUC@y.values[[1]] 
    return(c(AUC = AUC@y.values[[1]], AUC_replF = AUC_replF@y.values[[1]], AUC_diff = AUC_diff))
  })
  colnames(res_sep) = c(1:n)
  
  
  # Add col with permuted AUC vs original AUC
  res_sep = cbind(Permuted = c(ROC_PR_glm_concat_sig[[tissue_model]]$auc@y.values[[1]], AUC_permuted, 
                    AUC_permuted-ROC_PR_glm_concat_sig[[tissue_model]]$auc@y.values[[1]]), res_sep)
  
  
  
  
  ### All in one plot ###
  # ROC and PR plot concat
  print("plotting...")
  png(paste0("fig/featureExchange/ROC_PR_",tissue_model,"With", tissue_feature,"_",f,".png"), 
      height = 2000,width = 2000,pointsize=38)
  par(mfrow=c(2,2), cex.lab = 1.4, mar =  c(5, 4.2, 1.5, 0.5))

  # Scatterplot feature values between the selected two tissues
  plot(cbind(old_feature, repl_feature), xlab = paste0(p2P[f], " from ", t2T[tissue_model]), ylab = paste0(p2P[f], " from ", t2T[tissue_feature]), 
          pch = 16, cex = 0.4, col = alpha("black", 0.2), las = 1)
  lines(x = c(0,1), y = c(0,1), col = "gray35", lty = 2, lwd = 2)
  
  # Barplot AUC of segments 
  #barplot(res_sep[3,], ylab="AUC difference",xlab = "Decreasing feature difference", col = c("gray35", brewer.pal(n = 6, name = 'RdBu'))) 
  barplot(res_sep[3,], ylab="AUROC difference",xlab = "Decreasing feature difference", col = c("gray35", rep("gray75", 6)), las = 1) 
  

  # Boxplots of prediction difference per segment 
  boxplot(combPredsConcat$PredDiff ~ as.numeric(combPredsConcat$seg), 
          xlab = "Decreasing feature difference", ylab = "Prediction difference", col = rep("gray75", 6), las = 1)      #, col=brewer.pal(n = 6, name = 'RdBu')
  abline(h=0, col = "gray35", lty = 2, lwd = 2)
  
  # Barplot AUC, AUC_exchanged, AUC_permuted
  bar_auc = c(ROC_PR_glm_concat_sig[[tissue_model]]$auc@y.values[[1]], AUC_replF, AUC_permuted)
  barplot(bar_auc, ylab="AUROC", col = c(brewer.pal(n = 6, name = 'RdBu')[c(6,1)], "gray35"), 
          names =  c("Original", "Replaced", "Permuted"), las =1)
  dev.off()
})