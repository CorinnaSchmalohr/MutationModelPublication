tissue = "Skin"
load(paste0("temp/tempRes_", tissue, ".RData")) # res
load(paste0("data/MutTables/SomamutDB/", 
            tissue, "_WGS.RData")) # meta, Muts
res = cbind(Muts, prediction = res$pred);rm(Muts)
rownames(res) = paste0(res$chr, "_", res$pos)

# split Muts into 4 possible source papers and recompute performance
# match meta with predictions, using chr_pos label
res$paper = meta[rownames(res),"Paper"]

# pair each TP with a TN with the same sequence context
# or, in other words, assign the TNs with papers, sampled from the from TPs with the same context
load("data/processedData/pentamers.RData")
fivemerTranslator  = setNames(nm=c(fivemers[,1], fivemers[,2]),
                              object=c(rownames(fivemers), rownames(fivemers)))
for(context in rownames(fivemers)){
  res$paper[(res$context %in% fivemers[context,]) & res$mutated == 0] = 
    sample(res$paper[(res$context %in% fivemers[context,]) & res$mutated == 1])
}
table(res$paper)

# compute performance for each paper separately
library(ROCR)
perfByPaper = sapply(unique(res$paper), function(pap){
  subRes = res[res$paper == pap,]
  temp = prediction(pred = subRes$pred,labels = subRes$mutated)
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)

# visualize
# ROC
png("fig/healthyTissuesWGS/skin_perfByProject.png")
plot(NA, xlim = c(0,1), ylim = c(0,1),
     xlab = "FPR", ylab = "TPR", las = 1, main = tissue)
abline(0,1, lty = 2, col = "grey")
dumpVar = sapply(1:length(perfByPaper), function(i){
  lines(perfByPaper[[i]]$roc@x.values[[1]], perfByPaper[[i]]$roc@y.values[[1]], col = i+1)
})
legend("topleft", legend = names(perfByPaper), col = (1:length(perfByPaper)+1), lty = 1)
dev.off()
# the genomic landscapes of individual melanocytes... : paper says they also did WES. Check that This is reflected in annotation from somamutDB
# effects of psoriasis, ...:  is exome sequencing

