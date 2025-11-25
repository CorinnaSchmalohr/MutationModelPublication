library(Biostrings)
library(reshape2)
library(ggplot2)
library(ggforce)
library(ROCR)
library(gplots)
source("scripts/05_analysis/00_NamesAndColors.R")
load("data/processedData/dataInfos.RData")
plotEnding = "_20230905"


# prepare sequence contexts #####
bases = c("A", "C", "G", "T")
trimers = sort(apply(expand.grid(bases, bases, bases),
                      1,paste0, collapse = ""))
trimers = cbind(trimers,
                as.character(reverseComplement(DNAStringSet(trimers))))
trimers = t(apply(trimers,1,function(x){
  temp = which(substr(x,2,2) %in% c("C", "T"))
  return(c(x[temp], x[-temp]))
}))
trimers = unique(trimers)
rownames(trimers) = paste(trimers[,1], trimers[,2], sep = "/")
trimers = trimers[order(substr(trimers[,1], 2,2),  
                          substr(trimers[,1], 1,1), 
                          substr(trimers[,1], 3,3)),]
trim2context = setNames(nm=c(trimers[,1], trimers[,2]),
                        object=c(rownames(trimers), rownames(trimers)))
#####


# count how often each trimer appears and visualize ####
trimerFreqs = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  sub = substr(data$muts$context,start = 2,stop = 4)
  fr = table(factor(sub, levels = trimers[,1]))
  fr2 = table(factor(sub, levels = trimers[,2]))
  res = fr+fr2
  names(res) = rownames(trimers)
  return(res)
})
png(paste0("fig/ContextProportionsPerTissue.png"), width = 1600, height = 700, pointsize = 20)
trimerProp =apply(trimerFreqs, 2,function(x){x/sum(x)})
layout(rbind(1:10))
par(mar = c(0,1,2,0), oma = c(3.3,6,1,0))
dumpVar = sapply(tissues, function(tissue){
  barplot(rev(trimerProp[,tissue]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
          names.arg="",cex.names=0.5, cex.axis = 1.1,
          col = tissueCols[tissue],   main = t2T[tissue])
  if(tissue == tissues[1]){
    names = rev(rownames(trimerProp))
    axis(2, at = 1:length(names)-0.5, labels = names, las = 1, family = "mono", cex = 0.7)
  }
  abline(h=0, lty = 1, lwd = 2)
})
mtext("Proportion", side = 1, outer = T, line = 2)
mtext("Context", side = 2, outer = T, line = 4)
dev.off()
#####


# get prediction performance per trimer #####
load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")

perfTrimers = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  preds = predPerTissueRF[[tissue]]
  temp = do.call(rbind,preds)
  if( !identical(data$muts$mutated, (as.numeric(temp$label)-1))){
    stop("Predictions and mutations don't match")
  }
  mutationCalls = cbind(trimer = substr(data$muts[,"context"],2,4), temp)
  mutationCalls$trimer = trim2context[mutationCalls$trimer]
  mutationCallsPerTrimer = split(mutationCalls, mutationCalls$trimer)
  perfs = sapply(mutationCallsPerTrimer, function(x){
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  }, simplify = F)
  return(perfs)
}, simplify = F)
#####


# visualize performance #####
# ROC
png(paste0("fig/ContextWisePerformanceROC", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(0,1)
  dumpVar2 = sapply(perfTrimers[[tissue]], function(x){
    lines(x$roc@x.values[[1]], 
         x$roc@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "TPR", side = 2, line = 2)
  }
})
mtext(text = "FPR", side = 1, line = 0, outer = T)
dev.off()
# PR
png(paste0("fig/ContextWisePerformancePR", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(h = 0.5)
  dumpVar2 = sapply(perfTrimers[[tissue]], function(x){
    lines(x$pr@x.values[[1]], 
          x$pr@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "Precision", side = 2, line = 2)
  }
})
mtext(text = "Recall", side = 1, line = 0, outer = T)
dev.off()
# AUROC
AUCtrimer = sapply(perfTrimers, function(x){
  sapply(x, function(y){y$auc})
})
png(paste0("fig/ContextWisePerformanceAUROC", plotEnding, ".png"), 
    width = 800, height = 800, pointsize = 15)
heatmap(AUCtrimer, Rowv = NA, Colv = NA, labCol = t2T[tissues], 
        col = hcl.colors(100, palette = "viridis"),scale = "none")
dev.off()
# AUROC versus number of available datapoints
png(paste0("fig/ContextWisePerformanceAUROCvsNmuts", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(2,5), mar = c(3,3,2,1), oma = c(2,2,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(trimerFreqs[,tissue], AUCtrimer[,tissue], main=t2T[tissue], las = 1)
})
mtext("Number of positions", side = 1, outer = T)
mtext("AUROC", side = 2, outer = T)
dev.off()

png(paste0("fig/ContextWisePerformanceAUROCvsNmutsCombined", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(1,2))
plot(t(trimerFreqs), t(AUCtrimer),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Count", ylab = "AUROC")
legend("topright", col = tissueCols, legend = names(tissueCols),
       pch = 21, pt.bg = paste0(tissueCols, "65"), ncol = 2, cex = 0.95)
plot(t(trimerProp), t(AUCtrimer),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Proportion", ylab = "AUROC")
dev.off()
#####


# prepare mutation type #####
mutClass = rbind(c("C>A", "G>T"),
                 c("C>G", "G>C"),
                 c("C>T", "G>A"),
                 c("T>A", "A>T"),
                 c("T>C", "A>G"),
                 c("T>G", "A>C"))
rownames(mutClass) = mutClass[,1]
mutClassTranslator = setNames(nm=c(mutClass[,1], mutClass[,2]),
                              object=c(rownames(mutClass), rownames(mutClass)))
#####


# count how often each mutation class appears and visualize #####
typeFreqs = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  tempData = data$muts[data$muts$mutated == 1,] # get only TPs
  muttype = paste0(tempData$ref, ">", tempData$alt)
  muttype = mutClassTranslator[muttype]
  res = table(factor(muttype, levels = rownames(mutClass)))
  return(res)
})
typeProps =apply(typeFreqs, 2,function(x){x/sum(x)})

pdfAndPng(file = paste0("fig/MuttypeProportionsPerTissue"), 
          width = 16, height = 3, pngArgs = list(pointsize = 20), pdfArgs = list(pointsize=15),
          expr = expression({
            layout(rbind(1:10))
            par(mar = c(0,1,2,0), oma = c(3.3,6,1,1))
            dumpVar = sapply(tissues, function(tissue){
              barplot(rev(typeProps[,tissue]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
                      names.arg="",cex.names=0.5, cex.axis = 1.1,
                      col = tissueCols[tissue],  main = t2T[tissue])
              if(tissue == tissues[1]){
                names = rev(rownames(typeProps))
                axis(2, at = 1:length(names)-0.5, labels = names, las = 1, family = "mono", cex = 0.7)
              }
              abline(h=0, lty = 1, lwd = 2)
            })
            mtext("Proportion", side = 1, outer = T, line = 2)
            mtext("Mutation type", side = 2, outer = T, line = 4)
          }))
######


# get prediction performance per mutation Type #####
load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
load("data/processedData/pentamers.RData")
fivemerTranslator  = setNames(nm=c(fivemers[,1], fivemers[,2]),
                              object=c(rownames(fivemers), rownames(fivemers)))
perfMuttype = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  preds = predPerTissueRF[[tissue]]
  temp = do.call(rbind,preds)
  if( !identical(data$muts$mutated, (as.numeric(temp$label)-1))){
    stop("Predictions and mutations don't match")
  }
  mutationCalls = cbind(data$muts[,c("ref", "alt", "context")], temp)
  mutationCalls$mutType = mutClassTranslator[paste0(mutationCalls$ref, ">", mutationCalls$alt)]
  # pair each TP with a TN with the same sequence context
  # or, in other words, assign the TNs with mutation types from TPs with the same context
  for(context in rownames(fivemers)){
    mutationCalls$mutType[(mutationCalls$context %in% fivemers[context,]) & mutationCalls$label == 0] = 
      sample(mutationCalls$mutType[(mutationCalls$context %in% fivemers[context,]) & mutationCalls$label == 1])
  }
  
  mutationCallsPerMuttype = split(mutationCalls, mutationCalls$mutType)
  perfs = sapply(mutationCallsPerMuttype, function(x){
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  }, simplify = F)
  return(perfs)
}, simplify = F)
#####


# visualize prediction performance per mutation Type ####
# ROC
png(paste0("fig/MuttypeWisePerformanceROC", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(0,1)
  dumpVar2 = sapply(perfMuttype[[tissue]], function(x){
    lines(x$roc@x.values[[1]], 
          x$roc@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "TPR", side = 2, line = 2)
  }
})
mtext(text = "FPR", side = 1, line = 0, outer = T)
dev.off()
# PR
png(paste0("fig/MuttypeWisePerformancePR", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(h = 0.5)
  dumpVar2 = sapply(perfMuttype[[tissue]], function(x){
    lines(x$pr@x.values[[1]], 
          x$pr@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "Precision", side = 2, line = 2)
  }
})
mtext(text = "Recall", side = 1, line = 0, outer = T)
dev.off()
# AUROC
AUCmuttype = sapply(perfMuttype, function(x){
  sapply(x, function(y){y$auc})
})
pdfAndPng(file = paste0("fig/MuttypeWisePerformanceAUROC", plotEnding), 
          width = 8, height = 8, pngArgs = list(pointsize = 15), pdfArgs = list(pointsize=13), 
          expr = expression({
            heatmap.2(AUCmuttype, Rowv = NA, Colv = NA, dendrogram = "none",
                      labCol = t2T[tissues], trace = "none",scale = "none",
                      lhei = c(0.1, 4),lwid = c(0.1, 4),
                      cexRow = 2, cexCol = 2,margins = c(8,6),
                      srtCol = 45,adjCol = c(1,0.5),offsetCol = 0,key = F,
                      col = hcl.colors(100, palette = "viridis"),seq(from=0.5, to=0.7, length.out = 101)) #set the colors between 0.5 and 0.7, same as other plots
          }))
# 0.5-0.7
# AUROC versus number of available datapoints
png(paste0("fig/MuttypeWisePerformanceAUROCvsNmuts", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(2,5), mar = c(3,3,2,1), oma = c(2,2,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(typeFreqs[,tissue], AUCmuttype[,tissue], main=t2T[tissue], las = 1)
})
mtext("Number of positions", side = 1, outer = T)
mtext("AUROC", side = 2, outer = T)
dev.off()

png(paste0("fig/MuttypeWisePerformanceAUROCvsNmutsCombined", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(1,2))
plot(t(typeFreqs), t(AUCmuttype),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Count", ylab = "AUROC")
legend("topright", col = tissueCols, legend = names(tissueCols),
       pch = 21, pt.bg = paste0(tissueCols, "65"), ncol = 2, cex = 0.95)
plot(t(typeProps), t(AUCmuttype),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Proportion", ylab = "AUROC")
toMark = which(typeProps>0.4 | AUCmuttype>0.65, arr.ind = T)
text(typeProps[toMark],AUCmuttype[toMark],  
     labels=paste0("\n",rownames(toMark), "\n",t2T[colnames(AUCmuttype)[toMark[,2]]]),
     cex= 0.7, pos=2)
dev.off()
pdfAndPng(paste0("fig/MuttypeWisePerformanceAUROCvsPropMuts", plotEnding), 
    width = 10, height = 10, pngArgs = list(pointsize = 22), pdfArgs = list(pointsize=16), 
    expr = expression({
      plot(t(typeProps), t(AUCmuttype),  las = 1, col = tissueCols, pch = rep(c(0:6), each = length(tissues)),#pch = 21, 
            xlab = "Proportion per tissue", ylab = "AUC",)
      toMark = which(typeProps>0.4 | AUCmuttype>0.65, arr.ind = T)
      rect(xleft = -1, xright = 0.4, ybottom = -1, ytop = 0.65, border = "grey", lty = 2 )
      # abline(h=0.65, v = 0.4, lty = 2, col = "grey")
      text(typeProps[toMark],AUCmuttype[toMark],  
           labels=paste0("\n",rownames(toMark), "\n",t2T[colnames(AUCmuttype)[toMark[,2]]]),
           cex= 1, pos=2)
      temp1 = legend(x = 0.45,y=0.58, col = tissueCols, legend = t2T[names(tissueCols)],
             pch = 1, pt.bg = paste0(tissueCols, "65"), ncol = 2, bty = "n", text.width = 0.08)
      temp2 = legend(x=0.45,y=0.535, pch = c(0:6), pt.bg = "grey",col = "black", bty = "n",
             legend = rownames(typeProps), ncol = 2, text.width = 0.08)
      rect(xleft = temp2$rect$left, xright = temp2$rect$left+temp2$rect$w, ybottom = temp2$rect$top-temp2$rect$h, ytop = temp1$rect$top)
    }))
#####


# combination: trimer + mutation type = signature #####
bases = c("A", "C", "G", "T")
signatures = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
                                c(rep("C", 48), rep("T", 48)), ">",
                                c(rep(c("A", "G", "T"), each =16),
                                  rep(c("A", "C", "G"), each =16)), "]",
                                rep(bases, 24))),
                  paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
                                c(rep("G", 48), rep("A", 48)), ">",
                                c(rep(c("T", "C", "A"), each =16),
                                  rep(c("T", "G", "C"), each =16)), "]",
                                rep(rev(bases), 24))))
rownames(signatures) = paste(signatures[,1])
signatureTranslator= setNames(nm=c(signatures[,1], signatures[,2]),
                              object=c(rownames(signatures), rownames(signatures)))
signatureCounts = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  muts = data$muts[data$muts$mutated == 1,]
  sign = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
                substr(muts$context,4,4))
  sign2 = signatureTranslator[sign]
  res = table(factor(sign2,levels = signatures[,1]))
  array(res, dimnames = list(names(res)))
})
signatureProps =apply(signatureCounts, 2,function(x){x/sum(x)})

#####



# and visualize ####
png(paste0("fig/SignatureProportionsPerTissue.png"), width = 1600, height = 700, pointsize = 20)
layout(rbind(1:10))
par(mar = c(0,1,1,0), oma = c(3.3,5,1,0))
dumpVar = sapply(tissues, function(tissue){
  barplot(rev(signatureProps[,tissue]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
          names.arg="",cex.names=0.5, cex.axis = 1.1, main = t2T[tissue], 
          col = rep(rev(basesCol), each=16), xlim = c(0,0.16), xaxp = c(0, 0.1, 1))
  if(tissue == tissues[1]){
    classes = signatures[rownames(signatureProps),1]
    names = rev(paste0(substr(classes,1,1), substr(classes,3,3), substr(classes,7,7)))
    names2 = rev(paste0(" ", substr(classes,3,3), " "))
    text(x = 0, y = 0.99:96-0.5,  labels=names,family = "mono",
         xpd = NA, adj=1.1,  cex = 0.5)
    text(x = 0, y = 0.99:96-0.5,  labels=names2,family = "mono",
         xpd = NA, adj=1.1, col = rep(rev(basesCol), each=16), font = 2, cex = 0.5)
    text(x = -0.055, y = 1:6*16-8,
         labels = c("T>G", "T>C", "T>A" ,"C>T", "C>G", "C>A"),
         xpd = NA, adj=1, col = rev(basesCol), font = 2)
    segments(y0=1:6*16-16, y1 = 1:6*16, x0=-0.045, lwd = 4, 
             col = rev(basesCol), xpd = NA, lend=1, cex = 1.1)
    title(ylab = "Mutation type", mgp = c(5,1,0), xpd = NA, cex.lab = 1.1)
  }
  abline(h=seq(16,80,length.out = 5), lty = 2, lwd = 1.5)
  abline(h=0, lty = 1, lwd = 2)
  #abline(h=96, lty = 1, lwd = 1.8)
})
title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
dev.off()
#####


# get prediction performance per signature #####
load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
load("data/processedData/pentamers.RData")
fivemerTranslator  = setNames(nm=c(fivemers[,1], fivemers[,2]),
                              object=c(rownames(fivemers), rownames(fivemers)))
perfSignature = sapply(tissues, function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  preds = predPerTissueRF[[tissue]]
  temp = do.call(rbind,preds)
  if( !identical(data$muts$mutated, (as.numeric(temp$label)-1))){
    stop("Predictions and mutations don't match")
  }
  mutationCalls = cbind(data$muts[,c("ref", "alt", "context")], temp)
  tempSignature = paste0(substr(mutationCalls$context,2,2), 
                         "[",mutationCalls$ref, ">", mutationCalls$alt, "]",
                         substr(mutationCalls$context,4,4))
  mutationCalls$signature = signatureTranslator[tempSignature]
  # pair each TP with a TN with the same sequence context
  # or, in other words, assign the TNs with mutation types from TPs with the same context
  for(context in rownames(fivemers)){
    mutationCalls$signature[(mutationCalls$context %in% fivemers[context,]) & mutationCalls$label == 0] = 
      sample(mutationCalls$signature[(mutationCalls$context %in% fivemers[context,]) & mutationCalls$label == 1])
  }
  mutationCalls$signature = factor(mutationCalls$signature, levels = rownames(signatures))
  mutationCallsPerSignature = split(mutationCalls, mutationCalls$signature)
  perfs = sapply(mutationCallsPerSignature, function(x){
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  }, simplify = F)
  return(perfs)
}, simplify = F)
#####


# visualize prediction performance per signature #####
# ROC
png(paste0("fig/SignatureWisePerformanceROC", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(0,1)
  dumpVar2 = sapply(perfSignature[[tissue]], function(x){
    lines(x$roc@x.values[[1]], 
          x$roc@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "TPR", side = 2, line = 2)
  }
})
mtext(text = "FPR", side = 1, line = 0, outer = T)
dev.off()
# PR
png(paste0("fig/SignatureWisePerformancePR", plotEnding, ".png"), 
    width = 1600, height = 300, pointsize = 15)
layout(rbind(1:10))
par(mar = c(3,0.5,1,0), oma = c(2,3,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1), #xlab = "FPR", ylab = "TPR", 
       main = t2T[tissue], las = 1, mgp = c(3,1,0), yaxt = "none")
  abline(h = 0.5)
  dumpVar2 = sapply(perfSignature[[tissue]], function(x){
    lines(x$pr@x.values[[1]], 
          x$pr@y.values[[1]], col = paste0(tissueCols[tissue],"64"))
  })
  if(tissue == tissues[1]){
    axis(2, las = 1)
    mtext(text = "Precision", side = 2, line = 2)
  }
})
mtext(text = "Recall", side = 1, line = 0, outer = T)
dev.off()
# AUROC
AUCsignature = sapply(perfSignature, function(x){
  sapply(x, function(y){y$auc})
})
png(paste0("fig/SignatureWisePerformanceAUROC", plotEnding, ".png"), 
    width = 800, height = 800, pointsize = 15)
heatmap.2(AUCsignature, Rowv = F, Colv = F, dendrogram = "none",
          scale = "none",
          col = hcl.colors(100, palette = "viridis"),
          trace = "none",
          labCol = t2T[tissues],cexRow = 0.6,
          density.info = "none",key = F,
          RowSideColors = rep(basesCol, each=16),
          srtCol = 60, offsetCol= 0, offsetRow = 0,
          lmat = rbind(c(5,4,0), c(3,2,1)),
          lwid=c(0.1,4,0.1), lhei = c(0.01,4))
dev.off()

####
# AUROC versus number of available datapoints
png(paste0("fig/SignatureWisePerformanceAUROCvsNmuts", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(2,5), mar = c(3,3,2,1), oma = c(2,2,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(signatureCounts[,tissue], AUCsignature[,tissue], main=t2T[tissue], 
       las = 1)
})
mtext("Number of positions", side = 1, outer = T)
mtext("AUROC", side = 2, outer = T)
dev.off()
png(paste0("fig/SignatureWisePerformanceAUROCvsNmuts_logScale", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(2,5), mar = c(3,3,2,1), oma = c(2,2,0,0))
dumpVar = sapply(tissues, function(tissue){
  plot(signatureCounts[,tissue], AUCsignature[,tissue], main=t2T[tissue], 
       las = 1, log = "x", xlim = c(10,max(signatureCounts[,tissue])))
})
mtext("Number of positions (log axis)", side = 1, outer = T)
mtext("AUROC", side = 2, outer = T)
dev.off()

png(paste0("fig/SignatureWisePerformanceAUROCvsNmutsCombined", plotEnding, ".png"), 
    width = 1200, height = 600, pointsize = 15)
par(mfrow = c(1,2))
plot(t(signatureCounts), t(AUCsignature),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Count", ylab = "AUROC")
legend("bottomright", col = tissueCols, legend = names(tissueCols),
       pch = 21, pt.bg = paste0(tissueCols, "65"), ncol = 2, cex = 0.95)
plot(t(signatureProps), t(AUCsignature),  las = 1, col = tissueCols, pch = 21, 
     bg = paste0(tissueCols, "65"), xlab = "Proportion", ylab = "AUROC")
toMark = which(signatureProps>0.1 & AUCsignature>0.6, arr.ind = T)
text(signatureProps[toMark],AUCsignature[toMark],  
     labels=paste0("\n",rownames(toMark), "\n",t2T[colnames(AUCsignature)[toMark[,2]]]),
     cex= 0.7, pos=1)
dev.off()
# remove AUROCs with too little data
png(paste0("fig/SignatureWisePerformanceAUROC_filtered", plotEnding, ".png"), 
    width = 800, height = 800, pointsize = 15)
AUCsignatureFiltered = AUCsignature
AUCsignatureFiltered[signatureCounts <= 100] = NA
heatmap.2(AUCsignatureFiltered, Rowv = F, Colv = F, dendrogram = "none",
          scale = "none",
          col = hcl.colors(100, palette = "viridis"), na.color = "grey",
          trace = "none",
          labCol = t2T[tissues],cexRow = 0.6,
          density.info = "none",key = F,
          colCol = tissueCols, RowSideColors = rep(basesCol, each=16),
          srtCol = 60, offsetCol= 0, offsetRow = 0,
          lmat = rbind(c(5,4,0), c(3,2,1)),
          lwid=c(0.1,4,0.1), lhei = c(0.01,4))
dev.off()

# now with another column with the averages
AUCsignatureWMean = cbind(AUCsignatureFiltered,
                          Mean = rowMeans(AUCsignature, na.rm = T),
                          Median = apply(AUCsignature,1,median, na.rm = T))
png(paste0("fig/SignatureWisePerformanceAUROC_withMean", plotEnding, ".png"), 
    width = 800, height = 800, pointsize = 15)
heatmap.2(AUCsignatureWMean, Rowv = F, Colv = F, dendrogram = "none",
          scale = "none",
          col = hcl.colors(100, palette = "viridis"),
          trace = "none",
          labCol = c(t2T[tissues], "Mean", "Median"),cexRow = 0.6,
          density.info = "none",key = F,
          colCol = c(tissueCols, "black", "black"), RowSideColors = rep(basesCol, each=16),
          srtCol = 60, offsetCol= 0, offsetRow = 0,
          lmat = rbind(c(5,4,0), c(3,2,1)),
          lwid=c(0.1,4,0.1), lhei = c(0.01,4))
dev.off()
#####

