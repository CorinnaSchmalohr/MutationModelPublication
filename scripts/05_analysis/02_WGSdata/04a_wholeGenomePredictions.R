library(ROCR)
library(GenomicRanges)
tissues =  c("brain", "breast", "esophagus", "kidney", "liver", "ovary" ,
             "prostate", "skin")
methods = c("RFpreds", "contextPreds", "mult", "combOdds")

# cancer performance
perf = sapply(tissues, function(tissue){
  print(tissue)
  # load exome-wide prediction table
  load(paste0("data/Modeling/WholeGenomeData/combinedPredictions/combinedPredictions_", 
              tissue, "_chr","22",".RData")) #data
  if(all(is.na(data$healthyMuts)))
    data$healthyMuts = NULL
  data = mcols(data)
  data = na.omit(data)
  # compute performance for cancer data (training performance)
  cat("cancer data: ")
  perfCancer = sapply(methods, function(method){
    cat(method, ' ')
    perf = prediction(data[,method], labels = data$cancerMuts)
    auc = performance(perf,"auc")
    return(auc@y.values[[1]])
  }, simplify = T) ; cat('\n')
  # compute test performance (healthy tissues)
  cat("healthy data: ")
  
  perfHealthy = sapply(methods, function(method){
    cat(method, ' ')
    if(tissue == "ovary")
      return(NA)
    perf = prediction(data[,method], labels = data$healthyMuts)
    auc = performance(perf,"auc")
    return(auc@y.values[[1]])
  }, simplify = T); cat('\n')
  save(perfCancer, perfHealthy, 
       file = paste0("data/Modeling/WholeGenomeData/perf_Allmethods", tissue, ".RData"))
  return(list(perfCancer = perfCancer, perfHealthy = perfHealthy))
}, simplify = F)
save(perf, file = paste0("data/Modeling/WholeGenomeData/perf_Allmethods_", "allTissuesTogether", ".RData"))

# visualize:comparing RF, context, and combination for cancer and healthy mutations each
load(paste0("data/Modeling/WholeGenomeData/perf_Allmethods_", "allTissuesTogether", ".RData"))
dir.create("fig/wholeGenomePredictions", showWarnings = F)
source("scripts/05_analysis/00_NamesAndColors.R")
methodCols = setNames(rainbow(length(methods)), methods)



png("fig/wholeGenomePredictions/AUC_combinedMethods.png", 
    width = 800, height = 1200, pointsize = 30)
par(mfrow = c(4,2), mar = c(2,3,1,1), oma = c(1,1,1,10))
dumpVar = sapply(WGStissues, function(tissue){
  dat2plot = do.call(rbind,perf[[tissue]])
  temp = barplot(t(dat2plot), beside = T, las = 1, ylab = "AUC", 
                 col = methodCols, main = t2T[tissue], axisnames = F,
                 mgp = c(2,1,0)) #, ylim = c(0.35, max(aucs))
  axis(1,at = colMeans(temp), labels = c("Cancer data", "Healthy tissue data"), 
       line = -0.5, tick = F)
  if(tissue == WGStissues[2]){
    legend(x=11, y=0.7,fill = methodCols, legend = names(methodCols), xpd = NA)
  }
})
dev.off()
