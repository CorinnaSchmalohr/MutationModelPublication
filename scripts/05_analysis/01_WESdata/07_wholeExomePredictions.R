library(ROCR)
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung", "ovary", 
            "prostate", "skin") # "allTissues"
methods = c("RFpreds", "contextPreds", "mult", "mean", "combOdds", "LMcombination")

# cancer performance
perf = sapply(tissues, function(tissue){
  print(tissue)
  # load exome-wide prediction table
  load(paste0("data/Modeling/WholeExomeData/combinedPredictions/combinedPredictions_", tissue, ".RData")) #data
  if(all(is.na(data$healthyMuts)))
    data$healthyMuts = NULL
  data = na.omit(data)
  # compute performance for cancer data (training performance)
  cat("cancer data: ")
  perfCancer = sapply(methods, function(method){
    cat(method, ' ')
    perf = prediction(data[,method], labels = data$cancerMuts)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")
    return(list(roc = cbind(x = roc@x.values[[1]],
                            y = roc@y.values[[1]]),
                pr = cbind(x = pr@x.values[[1]],
                           y = pr@y.values[[1]]), 
                auc = auc@y.values[[1]]))
  }, simplify = F) ; cat('\n')
  # compute test performance (healthy tissues)
  cat("healthy data: ")
  if(all(is.na(data$healthyMuts))){
    perfHealthy = NA
  } else{
    perfHealthy = sapply(methods, function(method){
      cat(method, ' ')
      perf = prediction(data[,method], labels = data$healthyMuts)
      roc = performance(perf, "tpr", "fpr")
      pr = performance(perf,"prec", "rec")
      auc = performance(perf,"auc")
      return(list(roc = cbind(x = roc@x.values[[1]],
                              y = roc@y.values[[1]]),
                  pr = cbind(x = pr@x.values[[1]],
                             y = pr@y.values[[1]]), 
                  auc = auc@y.values[[1]]))
    }, simplify = F); cat('\n')
  }
  save(perfCancer, perfHealthy, 
       file = paste0("data/Modeling/WholeExomeData/perf_Allmethods", tissue, ".RData"))
  return(NULL)
})

# visualize: ROC and PR, comparing RF, context, and combination for cancer and healthy mutations each
dir.create("fig/wholeExomePredictions", showWarnings = F)
source("scripts/05_analysis/00_NamesAndColors.R")
methodCols = setNames(rainbow(length(methods)), methods)
png("fig/wholeExomePredictions/ROC_combinedMethods_cancer.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(5,2), mar = c(2,3,1,1), oma = c(1,1,1,10))
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  
  load(paste0("data/Modeling/WholeExomeData/perf_Allmethods", tissue, ".RData"))
  # cancer performance ROC
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.1,0.7,0), las = 1,xlab = "", ylab = "",main = t2T[tissue])
  abline(0,1,col = "grey", lty = 2)
  dumpVar2 = sapply(methods, function(method){
    subDat = perfCancer[[method]]$roc
    if(nrow(subDat) > 50000){
      subDat = subDat[sort(sample(nrow(subDat),50000)),]
    }
    lines(subDat, col = methodCols[method])
    return(NULL)
  })
  if(tissue == tissues[2]){
    legend(x=1.1,y=1, col = methodCols, lty = 1, legend = names(methodCols), xpd =NA)
  }
})
mtext()
dev.off()
png("fig/wholeExomePredictions/ROC_combinedMethods_healthy.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(5,2), mar = c(2,3,1,1))
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/Modeling/WholeExomeData/perf_Allmethods", tissue, ".RData"))
  # healthy performance ROC
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2.1,0.7,0), las = 1,xlab = "", ylab = "",main = t2T[tissue])
  abline(0,1,col = "grey", lty = 2)
  if(length(perfHealthy)>1){
    dumpVar2 = sapply(methods, function(method){
      subDat = perfHealthy[[method]]$roc
      if(nrow(subDat) > 100000){
        subDat = subDat[sort(sample(nrow(subDat),100000)),]
      }
      lines(subDat, col = methodCols[method])
      return(NULL)
    })
  } else{
    legend("topleft", col = methodCols, lty = 1, legend = names(methodCols))
  }
})
dev.off()



# collect AUC values
AUROCs = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/Modeling/WholeExomeData/perf_Allmethods", tissue, ".RData"))
  aucs = sapply(methods, function(method){
    res = c(perfCancer = perfCancer[[method]]$auc,
            perfHealthy = ifelse(all(is.na(perfHealthy)),NA,perfHealthy[[method]]$auc))
    return(res)

  })
  return(aucs)
}, simplify = F)


png("fig/wholeExomePredictions/AUC_combinedMethods.png", 
    width = 1200, height = 2000, pointsize = 30)
par(mfrow = c(5,2), mar = c(1,3,1,1), oma = c(1,1,1,10))
dumpVar = sapply(tissues, function(tissue){
  temp = barplot(t(AUROCS[[tissue]]), beside = T, las = 1, ylab = "AUC", 
                 col = methodCols, main = t2T[tissue], axisnames = F,
                 mgp = c(2,1,0), ylim = c(0.35, max(aucs)))
  axis(1,at = colMeans(temp), labels = c("Cancer data", "Healthy tissue data"), 
       line = -0.5, tick = F)
  if(tissue == tissues[2]){
    legend(x=15, y=0.7,fill = methodCols, legend = names(methodCols), xpd = NA)
  }
})
dev.off()
