.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ggplot2)
library(data.table)
library(ROCR)
library(ranger)
library(berryFunctions) # for smallPlot
library(corrplot)
library(sinaplot)
library(RColorBrewer)
library(plotrix) # for point labels
library(readxl)
source("lib/general_function.R")
source("scripts/05_analysis/00_NamesAndColors.R")
methodCols = methodCols[1:2]
tissues = c("brain", "breast", "esophagus", "kidney", "liver", "ovary", "prostate", "skin")
t2T = t2T[tissues]
dir.create("fig/modelEvaluationWGS", showWarnings=F)
plotEnding = "_20240902"


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
p2PWGS = setNames(predictorOrder$Name, predictorOrder$abbreviation)
predictorGroups = split(predictorOrder$abbreviation, predictorOrder$Group)
rm(tab, predictorOrder)
#####


# load data #####
load("data/processedData/WholeGenomedataInfos.RData")
load("data/Modeling/WholeGenomeData/RF/predPerTissueRF.RData")
load("data/Modeling/WholeGenomeData/GLM/predPerTissueGLM_sig.RData")
# load("data/Modeling/WholeGenomeData/Lasso/predPerTissueLasso.RData")
load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_perChr.RData")
load("data/Modeling/WholeGenomeData/GLM/ROC_PR_glm_perChr_sig.RData")
# load("data/Modeling/WholeGenomeData/Lasso/ROC_PR_lasso_perChr.RData")
load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_concat.RData")
load("data/Modeling/WholeGenomeData/GLM/ROC_PR_glm_concat_sig.RData")
# load("data/Modeling/WholeGenomeData/Lasso/ROC_PR_lasso_concat.RData")
load("data/Modeling/WholeGenomeData/RF/RF_imps.RData")
load("data/Modeling/WholeGenomeData/GLM/GLM_impsAndPvals_sig.RData")
# load("data/Modeling/WholeGenomeData/Lasso/Lasso_impAndStab.RData")
load("data/MutTables/WholeGenomeData/nPerTissue.RData")
load("data/Modeling/WholeGenomeData/RF/ROC_PR_RF_train.RData")
load("data/Modeling/WholeGenomeData/RF/ROC_PR_glm_sig_train.RData")
# load("data/Modeling/WholeGenomeData/RF/ROC_PR_lasso_sig_train.RData")
#####



# # Corrplots ######
# load("data/processedData/WholeGenomedataInfos.RData")
# dumpVar = sapply(tissues, function(tissue){
#   pdf(paste0("fig/modelEvaluationWGS/corrplot_", tissue, plotEnding, ".pdf"),
#       width=1550, height=1300, pointsize = 25)
#   png(paste0("fig/modelEvaluationWGS/corrplot_", tissue, plotEnding, ".png"),
#       width=1550, height=1300, pointsize = 25)
#   cors = dataInfos[[tissue]]$cors
#   rownames(cors) = p2PWGS[rownames(cors)]
#   colnames(cors) = p2PWGS[colnames(cors)]
#   corrplot(corr = cors, tl.cex=0.4, tl.col="black")
#   dev.off()
#   dev.off()
# })
# #####
# 
# # Mutation overview #####
# print("mutation overview")
# # prepare for upper part: mutation counts per sample, ordered
# sortByMedian = nPerTissue[order(sapply(nPerTissue,median))]
# # prepare for lower part:
# # counts of mutation types
# mutClass = list(c("T>A", "A>T"),
#                 c("T>C", "A>G"),
#                 c("T>G", "A>C"),
#                 c("C>A", "G>T"),
#                 c("C>G", "G>C"),
#                 c("C>T", "G>A"))
# names(mutClass) = sapply(mutClass, paste, collapse = "/")
# mutTypes = sapply(names(sortByMedian), function(tissue){
#   load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
#               tissue, "_mapped.RData"))
#   muts = data$muts[data$muts$mutated == 1,]
#   temp = paste(muts$ref, muts$alt, sep=">")
#   sapply(mutClass, function(x){
#     sum(temp %in% x)
#   })
# })
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
# trimTypes = sapply(names(sortByMedian), function(tissue){
#   load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
#               tissue, "_mapped.RData"))
#   muts = data$muts[data$muts$mutated == 1,]
#   temp = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
#                 substr(muts$context,4,4))
#   apply(trimClass,1, function(x){
#     sum(temp %in% x)
#   })
# })
# save(trimClass, trimTypes, file = "data/MutTables/WholeGenomeData/WGSMuts_trimClass_trimTypes.RData")
# # load("data/MutTables/WholeGenomeData/WGSMuts_trimClass_trimTypes.RData")
# 
# 
# 
# 
# # plot
# print("plotting")
# png(paste0("fig/MutationOverview_WGS", plotEnding, ".png"),
#     res=200, pointsize=18, width = 2900, height = 1500, )
# layout(rbind(c(1,1,1,1,1,1,1,1),
#              c(2,3,4,5,6,7,8,9),
#              c(2,3,4,5,6,7,8,9),
#              c(2,3,4,5,6,7,8,9)))
# # par(mar = c(0.1,0.15,0,0), oma = c(3.3,5,2,0))
# par(mar = c(0.1,1,0,0), oma = c(3.3,5,2,0))
# # upper part: mutation counts per sample, ordered
# sortByMedian = nPerTissue[order(sapply(nPerTissue,median))]
# plot(NA, xlim = c(1,16), ylim = c(1,max(unlist(nPerTissue))), log = "y",bty="l",
#      ylab = "Mutations per tumor", xaxt = "n", las = 1, xlab = "", xpd = NA,
#      mgp = c(5,1,0), cex.lab = 1.1,cex.axis = 1.1)
# dumpVar = sapply(1:length(sortByMedian), function(i){
#   tissue= names(sortByMedian)[i]
#   x = sortByMedian[[i]]
#   x = sort(x)
#   pos = seq(0,1,length.out = length(x))
#   points(i*2-1+pos,x, col = tissueCols[tissue], pch = 1, cex=0.7)
#   segments(x0=i*2-1,x1=i*2, y0 = median(x), lwd = 2, lty = 2)
#   # return(cbind(i*2-1+pos,x))
# })
# axis(3,at =1:length(sortByMedian)*2-0.5, labels= t2T[names(sortByMedian)], tick=F,
#      mgp = c(2,0.5,0), cex.axis  = 1.1)
# # lower part: frequencies of mutation types
# par(mar = c(0,1,0,0))
# temp =apply(trimTypes, 2,function(x){x/sum(x)})
# dumpVar = sapply(1:ncol(temp), function(i){
#   barplot(rev(temp[,i]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
#           names.arg="",cex.names=0.5, cex.axis = 1.1,
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
# title(xlab = "Proportion of variants", outer = T, mgp = c(2,1,0), cex.lab = 1.1)
# dev.off()
# #####
# 
# # Feature effect directions ######
# preds = names(p2PWGS)
# for(tissue in tissues){
#   print(tissue)
#   load(paste0("data/MutTables/WholeGenomeData/WGSMuts_",
#               tissue, "_mapped_processed.RData"))
#   pdf(paste0("fig/predictor_effectDirection_WGS_", tissue, plotEnding, ".pdf"),
#       width = 20,height = 15,pointsize=20)
#   par(mfrow = c(5,8), mar = c(2.5,3,3,0.5))
#   sapply(preds, function(x){
#     title = strwrap(p2PWGS[x], width = 18)
#     if(! x %in% colnames(dat)){
#       # plot(0,type='n',axes=FALSE, main = title, xlab = "", ylab = "", cex.main=0.8)
#     }else{
#       nLevels = length(unique(dat[,x]))
#       if(nLevels == 2){
#         pd = cbind(table(dat[,x][dat$mutated == 0])/sum(dat$mutated == 0), NA)
#         barplot(pd, names.arg=c("TN", "TP"),
#                 las = 1, main = title, cex.main=0.8,
#                 col = c(rgb(238,44,44,100, maxColorValue = 255),
#                         rgb(238,44,44, maxColorValue = 255)), legend.text =F)
#         pd = cbind(NA, table(dat[,x][dat$mutated == 1])/sum(dat$mutated == 1))
#         barplot(pd, col = c(rgb(0,0,139,100, maxColorValue = 255),
#                             rgb(0,0,139, maxColorValue = 255)),
#                 add= T, yaxt = "n", width=c(1,1))
#       } else if (nLevels <=5){
#         barplot(cbind(table(dat[,x][dat$mutated == 0]),
#                       table(dat[,x][dat$mutated == 1]), NA),
#                 names.arg=c("TN", "TP", ""),
#                 las = 1, main = title, cex.main=0.8,
#                 col = rainbow(nLevels), legend.text=T,
#                 args.legend = list(bty = "n"))
#       }
#       else if(is.numeric(dat[,x])){
#         temp = density(dat[,x][dat$mutated == 1])
#         temp2 = density(dat[,x][dat$mutated == 0])
#         plot(temp, col = rgb(238,44,44, maxColorValue = 255), lwd = 2.5,
#              main = title, xlab = "", las = 1,cex.main=0.8,
#              ylim = c(0,max(c(temp$y, temp2$y))))
#         polygon(temp, col = rgb(238,44,44,alpha = 75, maxColorValue = 255), border = NA)
#         lines(temp2, col = rgb(0,0,139, maxColorValue = 255), lwd = 2.5, lty = 2)
#         polygon(temp2, col = rgb(0,0,139,alpha = 75, maxColorValue = 255), border = NA)
#       }  else
#         stop("something went wrong")
#     }
#   })
#   dev.off()
# }
# #####


# # Comparison of RF, GLM and LASSO detailed plot #####
# print("Comparison of RF, GLM a detailed plot")
# 
# png(paste0("fig/modelEvaluationWGS/compareMethodperformance", plotEnding,".png"),
#     height=1600, width=1300, pointsize=30)#, res = 200)
# par(mfrow = c(length(tissues),4), mar = c(0.5,4,0.1,0.1),
#     mgp = c(2.5,1,0), oma = c(3,1,0.05,0.05))
# plotDump = sapply(tissues, function(tissue){
#   print(tissue)
#   # roc
#   print("roc")
#   plot(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
#        ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],lwd = 3,
#        col = methodCols["RF"], las = 1,xlab = "",
#        ylab = "", type = "l", xaxt = "n",yaxt = "n")
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = 1, padj = 0.7, cex.axis = 1.4)
#   abline(0,1, col = "grey", lty = 2)
#   lines(x = ROC_PR_glm_concat_sig[[tissue]]$roc@x.values[[1]],
#         y = ROC_PR_glm_concat_sig[[tissue]]$roc@y.values[[1]],
#         col = methodCols["GLM"], lty = 2, lwd = 3)
#   mtext(side=2, text=t2T[tissue], line = 3.7)
#   # add roc axis labels
#   print("roc axis")
#   if(tissue == tail(tissues,1)){
#     axis(1, mgp = c(2,0.7,0), las = 1,
#          at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis=1.4)
#     axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis=1.4)
#     title(xlab = "FPR", line = 2, xpd = NA, cex.lab = 1.4)
#   }
#   if(tissue == tissues[ceiling(length(tissues)/2)]){
#     title(ylab = "TPR", xpd = NA, mgp = c(2.5,1,0), cex.lab = 1.4)
#   }
#   # auroc
#   print("auroc")
#   AUCs = list(RF = sapply(ROC_PR_RF_perChr[[tissue]],function(x){x$auc}),
#               GLM = sapply(ROC_PR_glm_perChr_sig[[tissue]],function(x){x$auc}))
#               # SL = sapply(ROC_PR_lasso_perChr[[tissue]],function(x){x$auc}))
#   sinaplot(AUCs, adjust =0.1, maxwidth = 0.5, col = methodCols,
#            xaxt = "n", yaxt = "n")
#   axis(2, mgp = c(2,0.7,0), las = 1,cex.axis = 1.4,
#        at = par("yaxp")[1:2], lwd = 0, lwd.ticks=1)
#   AUCs = c(RF = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
#            GLM = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]])
#   
#   points(AUCs, col = methodCols, pch = 19, cex = 2)
#   # add auroc axis labels
#   print("auroc axis")
#   if(tissue == tail(tissues,1)){
#     # mtext(text=names(methodCols), side=1,  at=1:3, cex.lab = 1.4)
#     axis(1,at = 1:2, labels=names(methodCols), mgp = c(2,0.7,0),
#          cex.axis=1.4)
#   }
#   if(tissue == tissues[ceiling(length(tissues)/2)]){
#     title(ylab = "AUC", xpd = NA, mgp = c(2.5,1,0), cex.lab = 1.4)
#   }
#   #pr
#   print("pr")
#   plot(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
#        ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], lwd = 3,
#        las = 1,xlab = "", xaxt = "n", type = "l", col= methodCols["RF"],
#        ylab = "",  ylim = c(0,1), mgp =  c(2,0.7,0),  xaxt = "n",yaxt = "n")
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = 1, padj = 0.7, cex.axis = 1.4)
#   lines(x = ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]],
#         y = ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]],
#         col = methodCols["GLM"], lty = 2, lwd = 3)
# 
#   rect(xleft = -0.005, xright=0.02, ybottom=0.5, ytop=1, lwd = 1.5)
#   # add pr axis labels
#   print("pr axis")
#   if(tissue == tail(tissues,1)){
#     axis(1, mgp = c(2,0.7,0), las = 1,
#          at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
#     axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis = 1.4)
#     title(xlab = "Recall", line = 2, xpd = NA, cex.lab = 1.4)
#   }
#   if(tissue == tissues[ceiling(length(tissues)/2)]){
#     title(ylab = "Precision", xpd = NA,
#           mgp =  c(2.5,0.7,0), cex.lab = 1.4)
#   }
#   # create inset with zoom
#   print("inset")
#   lim = par("plt")
#   xlims = lim[2]-lim[1]
#   ylims = lim[4]-lim[3]
#   inset = 0.05
#   smallPlot(expr={
#     plot(x = ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
#          y = ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
#          xlab = "", ylab = "", xaxt = "n", yaxt = "n",
#          xlim = c(0,0.02),ylim = c(0.5,1),
#          lwd = 3, type = "l", col = methodCols["RF"])
#     lines(x = ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]],
#           y = ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]],
#           col = methodCols["GLM"], lty = 2, lwd = 3)
# 
#     box(lwd = 2)},
#     x1 = lim[1]+0.3*xlims, x2 = lim[2]-0.05*xlims,
#     y1 = lim[3]+0.05*ylims, y2 = lim[4]-0.3*ylims, xpd = F,
#     mar = c(0,0,0,0), border = "transparent")
# 
#   # violin of predictions
#   print("violin")
#   rfpreds = do.call(rbind,predPerTissueRF[[tissue]])
#   glmpreds = do.call(rbind,predPerTissueGLMsig[[tissue]])
#   if(nrow(rfpreds)>1000000){
#     rfpreds = rfpreds[sample(1:nrow(rfpreds), 1000000),]
#     glmpreds = glmpreds[sample(1:nrow(glmpreds), 1000000),]
#   }
# 
#   rfpredssplit = split(rfpreds$pred, rfpreds$label)
#   glmpredssplit = split(glmpreds$pred, glmpreds$label)
#   plotDat = list("TN.RF" = rfpredssplit$`0`,
#                    "TP.RF" = rfpredssplit$`1`,
#                    "TN.GLM" = glmpredssplit$`0`,
#                    "TP.GLM" = glmpredssplit$`1`)
#   tempCols = RColorBrewer::brewer.pal(3,"Set2")
#   print("plotting")
#   sinaplot(plotDat, col = c(tempCols,methodCols)[c(1,4,2,5,3,6)],
#            xaxt = "n",yaxt = "n",las = 1,  ylab = "", cex = 0.8,
#            ylim= c(0,1))
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
#   axis(2, mgp = c(2,0.7,0), las = 1,
#        at = 1, padj = 0.7, cex.axis = 1.4)
#   abline(h=0.5, col = "grey", lty = 2, lwd = 2)
#   abline(v=c(2.5,4.5))
#   boxplot(plotDat, col = c("grey80", "grey40"),
#           add = T, boxwex = 0.2, outline = F, #col = rgb(0,0,0,alpha = 0),
#           ann = F, yaxt = "n", xaxt = "n", mgp = c(2,0.7,0))
#   text(x= c(1.5,3.5), y=0,labels=names(methodCols), cex = 1.2)
#   if(tissue == tail(tissues,1)){
#     axis(1,at = 1:4, labels=rep(c("0","1"),2), mgp = c(2,0.7,0),
#          cex.axis=1.2)
#     title(xlab = "True labels", line = 2, xpd = NA,
#           mgp = c(2,0.7,0), cex.lab = 1.4)
#   }
#   if(tissue == tissues[ceiling(length(tissues)/2)]){
#     title(ylab = "Prediction", xpd = NA,
#           mgp =  c(2.5,0.7,0), cex.lab = 1.4)
#   }
#   print("done")
# })
# dev.off()
# #####
# 
# 
# # # Comparison of AUC between RF, GLM and LASSO #####
# print("Comparison of AUC between RF, GLM and LASSO ")
# 
# AUCcollection = sapply(tissues, function(tissue){
#   c(RF = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
#     GLM = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]])
# })
# png(paste0("fig/modelEvaluationWGS/AUCbetweenMethods_allTissues",
#            plotEnding, ".png"), width = 1200, height = 500, pointsize = 20)
# par(mar = c(3,4,1,4))
# temp = barplot(AUCcollection-0.5, beside = T, las = 1, ylab = "AUC",
#                col = methodCols, axisnames = F, yaxt = "n", names.arg = t2T[tissues])
# axis(2, at = axTicks(2), labels = c(0,axTicks(2)[-1]+0.5), las = 1)
# axis.break(axis = 2, breakpos = 0.01, style = "gap", bgcol = "grey")
# title(ylab = "AUC")
# legend("topleft", legend = names(methodCols),
#        fill = methodCols, xpd = NA, bty = "n")
# mtext(t2T[tissues], side = 1, at = colMeans(temp))
# dev.off()
# #####
# # 
# # Train vs. test performance of RF, GLM and LASSO #####
# print("Train vs. test performance of RF, GLM and LASSO")
# RFaucs = sapply(tissues, function(tissue){
#   c(train = ROC_PR_RF_train[[tissue]]$auc@y.values[[1]],
#       test = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]])
# })
# GLMaucs = sapply(tissues, function(tissue){
#   c(train = ROC_PR_glm_sig_train[[tissue]]$auc@y.values[[1]],
#     test = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]])
# })
# png(paste0("fig/modelEvaluationWGS/TrainVsTestAUC",
#            plotEnding, ".png"), width = 1200, height = 500, pointsize = 25)
# par(mfrow = c(2,1), mar = c(1,4,0.5,0.5), oma = c(2,4,2,0))
# barplot(RFaucs, beside = T, las = 1,
#         ylab = "AUC", xaxt = "n",
#         col = rep(tissueCols, each = 2),
#         density = rep(c(NA, 20),ncol(RFaucs)),
#         legend.text = c("Training set", "Test set"),
#         args.legend = list(x="topright", density=c(NA,40), ncol = 2,
#                            inset=c(0.05,-0.43), xpd = NA, fill = "black"))
# mtext(text = "RF", side = 2, las = 1,  line = 5,) 
# abline(h=0)
# barplot(GLMaucs, beside = T, las = 1,
#         ylab = "AUC", names.arg = t2T[tissues],
#         col = rep(tissueCols, each = 2),
#         density = rep(c(NA, 20),ncol(GLMaucs)))
# mtext(text = "GLM", side = 2, las = 1,  line = 5) 
# abline(h=0)
# dev.off()
# #####
# 
# # ROC over chrs for all tissues #####
# print("ROC over chrs for all tissues")
# png(paste0("fig/modelEvaluationWGS/RF_ROCsOverChrsAllTissues",plotEnding,".png"),
#     width=1200, height=400, res=150)
# par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,4))
# plotDump = sapply(tissues, function(tissue){
#   plot(NA, xlim = c(0,1), ylim = c(0,1),
#        mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
#        xaxt = "n", yaxt = "n")
#   abline(0,1, col = "grey", lty = 2)
#   dump = sapply(names(chrCols), function(cr){
#     plot(ROC_PR_RF_perChr[[tissue]][[cr]]$roc,
#          col = rgb(0,0,0,0.5), add = T)
#   })
#   lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
#         ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],
#         col = tissueCols[tissue], lwd = 2.5)
#   text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0,0))
#   if(tissue %in% c("brain", "liver")){
#     axis(2, las = 1, mgp = c(2,0.7,0))
#   }
#   if(tissue %in% tissues[5:8]){
#     axis(1, mgp = c(2,0.7,0))
#   }
# })
# title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T)
# title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T)
# legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
#        col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
# dev.off()
# #####
# 
# 
# # PR over chrs for all tissues #####
# print("PR over chrs for all tissues")
# png(paste0("fig/modelEvaluationWGS/RF_PRsOverChrsAllTissues", plotEnding,".png"),
#     width=1200, height=400, res=150)
# par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,4))
# plotDump = sapply(tissues, function(tissue){
#   plot(NA, xlim = c(0,1), ylim = c(0,1),
#        mgp = c(2,0.8,0), las = 1,xlab = "", ylab = "",
#        xaxt = "n", yaxt = "n")
#   dump = sapply(names(chrCols), function(cr){
#     plot(ROC_PR_RF_perChr[[tissue]][[cr]]$pr,
#          col = rgb(0,0,0,0.5), add = T)
#   })
#   plot(ROC_PR_RF_concat[[tissue]]$pr,
#        col = tissueCols[tissue], add = T, lwd = 2.5, main = "")
#   text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0.5,0))
#   if(tissue %in% c("brain", "liver")){
#     axis(2, las = 1, mgp = c(2,0.7,0))
#   }
#   if(tissue %in% tissues[5:8]){
#     axis(1, mgp = c(2,0.7,0))
#   }
# })
# title(xlab = "Recall", mgp = c(1.8,0.7,0), xpd = NA, outer = T)
# title(ylab = "Precision", mgp = c(2,0.7,0), xpd = NA, outer = T)
# legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
#        col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
# dev.off()
# #####
# 
# # AUC versus chromosome features #####
# print("AUC versus chromosome features")
# load("data/processedData/chrLengths.RData")
# 
# png(paste0("fig/modelEvaluationWGS/RF_chromAUCsFeatures", plotEnding,".png"),
#     width=1200, height=400, pointsize=15)
# par(mfrow = c(2,4),mar = c(4,3,2,0.5), oma = c(1,2,0.5,0.5))
# # AUC versus training data size
# dumpVar = sapply(tissues, function(tissue){
#   AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
#   nMuts = dataInfos[[tissue]]$nMutsPerChr/1000
#   plot(nMuts, AUCs, las = 1, mgp = c(2,0.7,0),
#        ylab = "", xlab = "", main = t2T[tissue])
#   # mod = lm(AUCs~nMuts)
#   # abline(mod, col = "grey", lty = 2, lwd = 2)
# })
# # axis("AUC")
# mtext(text = "AUC", side = 2, line = 0, outer = T)
# mtext(text = "Chromsome data size (1000 positions)", side = 1, line = 0, outer = T)
# dev.off()
# #####
# 
# 
# # Overview for tissue comparison #####
# print("Overview for tissue comparison ")
# png(paste0("fig/modelEvaluationWGS/compareRFperformanceTissues", plotEnding, ".png"),
#     width=2400, height=800, pointsize=35)
# par(mfrow = c(1,3), mar = c(4,4,0.1,0.1), mgp = c(2.5,1,0))
# # roc
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      xlab = "False positive rate",
#      ylab = "True positive rate", las = 1)
# abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
#         ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],
#         col = tissueCols[tissue], lwd = 2)
# })
# legend("bottomright", col=tissueCols[names(ROC_PR_RF_concat)],
#        legend=t2T[names(ROC_PR_RF_concat)], lwd = 2)
# # pr
# plot(NA, xlim = c(0,1), ylim = c(0,1),
#      xlab = "Recall",
#      ylab = "Precision", las = 1)
# plotDump = sapply(tissues, function(tissue){
#   lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
#         ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
#         col = tissueCols[tissue], lwd = 2)
# })
# ## AUCs vs nMuts ##
AUCs = sapply(tissues, function(tissue){
  ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]]
})
AUCsPerChr = sapply(tissues, function(tissue){
  sapply(ROC_PR_RF_perChr[[tissue]], function(x){
    x$auc
  })
})
AUCsSD = apply(AUCsPerChr,2,sd)
nMuts = sapply(dataInfos, function(x){x$nMuts})
# 
# plot(nMuts, AUCs, las = 1,
#      bg = tissueCols[names(ROC_PR_RF_concat)], pch = 21, cex = 1.5,
#      ylim = c(min(colMeans(AUCsPerChr)-AUCsSD),
#               max(colMeans(AUCsPerChr)+AUCsSD)),
#      xlim = c(min(nMuts), max(nMuts)+50000),
#      xlab = "n positions", ylab = "AUROC")
# abline(lm(colMeans(AUCsPerChr) ~ nMuts), col = "grey", lty = 2)
# arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
#        y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
#        angle = 90, length = 0.1)
# points(nMuts, AUCs,
#        bg = tissueCols[names(ROC_PR_RF_concat)], pch = 21, cex = 1.5)
# # points(x=nMuts, y = AUCs, pch = 4, bg = tissueCols, cex = 1.5)
# text(x=nMuts, y = AUCs, labels=t2T[tissues], pos=1)
# dev.off()

pdfAndPng(file = paste0("fig/modelEvaluationWGS/compareRFAUCTissues", plotEnding),
          width=8, height=8, pngArgs = list(pointsize=20), 
          pdfArgs = list(pointsize=15), 
          expr = expression({
            plot(nMuts, AUCs, las = 1,
                 bg = tissueCols[names(ROC_PR_RF_concat)], pch = 21, cex = 1.5,
                 ylim = c(min(colMeans(AUCsPerChr)-AUCsSD),
                          max(colMeans(AUCsPerChr)+AUCsSD)),
                 xlim = c(min(nMuts), max(nMuts)+50000),
                 xlab = "n Positions", ylab = "AUC")
            abline(lm(colMeans(AUCsPerChr) ~ nMuts), col = "grey", lty = 2)
            arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
                   y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
                   angle = 90, length = 0.1)
            points(nMuts, AUCs,
                   bg = tissueCols[names(ROC_PR_RF_concat)], pch = 21, cex = 1.5)
            # points(x=nMuts, y = AUCs, pch = 4, bg = tissueCols, cex = 1.5)
            text(x=nMuts, y = AUCs, labels=t2T[tissues], pos=1)
          }))

# #####
# 
# 
# # TODO Mbp-wise performance #####
# # only makes sense for WGS
# #####
# 
# # Predictor importances over chromosomes RF ######
# print("Predictor importances over chromosomes RF")
# # rf gini importance
rf_imps = do.call(rbind,sapply(tissues, function(tissue){
  load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_",
              "finalModel.RData")) # rf, importance
  temp = cbind(rf_gini[[tissue]],final = importance[rownames(rf_gini[[tissue]])])
  temp =as.data.table(temp, keep.rownames = "feature")
  temp$feature = factor(temp$feature, levels = temp$feature)
  temp = melt(temp, id.vars = "feature",
              variable.name = "chromosome", value.name = "gini")
  temp$tissue = t2T[tissue]
  temp$gini_scaled = scale(temp$gini, center = F)
  rm(rf, importance); gc()
  return(temp)
}, simplify=F))
rf_imps$group = factor(p2GWGS[rf_imps$feature], levels = rev(unique(p2GWGS)))
save(rf_imps, file = "data/Modeling/WholeGenomeData/RF/rf_imps.RData")
# # load("data/Modeling/WholeGenomeData/RF/rf_imps.RData")
# ggplot(rf_imps, aes(x = chromosome, y = feature, fill = gini_scaled)) +
#   geom_raster() +
#   scale_fill_gradient(low="grey90", high="red",na.value="grey") +
#   labs(y = "Predictor", x = "Chromosome") +
#   labs(fill = "Scaled Gini\nImportance")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2PWGS)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
#         axis.title=element_text(size=14,face="bold")) +
#   facet_wrap(~tissue, nrow = 1)
# ggsave(paste0("fig/modelEvaluationWGS/RFgini_perChr", plotEnding, ".png"),
#        height=8, width=13)
# 
# ######
# 
# 
# # Predictor importances over chromosomes GLM post-selection ######
# print("Predictor importances over chromosomes GLM post-selection")
# # coefficients
# # collect coefficients, add info if significant
# glm_coeffs_sig= do.call(rbind,sapply(tissues, function(tissue){
#   load(paste0("data/Modeling/WholeGenomeData/GLM/", tissue, "_sig.RData"))
#   # coefficients
#   temp = cbind(glm_imps_sig[[tissue]],final = logR$coefficients[names(p2PWGS)])
#   temp = as.data.table(temp, keep.rownames = "feature")
#   temp$feature = factor(temp$feature, levels = temp$feature)
#   temp = melt(temp, id.vars = "feature",
#               variable.name = "chromosome", value.name = "coeff")
#   # and pvalues
#   temp2 = cbind(glm_pvals_sig[[tissue]],final = coef(summary(logR))[,4][-1][names(p2PWGS)])
#   temp2 = as.data.table(temp2, keep.rownames = "feature")
#   temp2$feature = factor(temp2$feature, levels = temp2$feature)
#   temp2 = melt(temp2, id.vars = "feature",
#                variable.name = "chromosome", value.name = "pval")
#   # add pvalues to coefficient table
#   temp$pval = temp2$pval
#   temp$tissue = t2T[tissue]
#   temp$coeff_scaled = scale(temp$coeff, center = F)
#   return(temp)
# },simplify=F))
# ggplot(glm_coeffs_sig, aes(x = chromosome, y = feature, fill = coeff)) +
#   geom_raster() +
#   scale_fill_gradient2(
#     low = "blue",mid = "gray90",high = "red",
#     breaks = c(min(glm_coeffs_sig$coeff, na.rm = T), 0,
#                max(glm_coeffs_sig$coeff, na.rm = T)),
#     labels = c(round(min(glm_coeffs_sig$coeff, na.rm = T), digits = -1), 0,
#                round(max(glm_coeffs_sig$coeff, na.rm = T), digits = -1))) +
#   labs(y = "Predictor", x = "Chromosome") +
#   labs(fill = "Coefficient")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2PWGS)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
#         axis.title=element_text(size=14,face="bold")) +
#   facet_wrap(c("tissue"), nrow = 1)
# ggsave(paste0("fig/modelEvaluationWGS/GLMcoeffSig_perChr", plotEnding, ".png"),
#        height=8, width=13)
# # pvalues
# glm_coeffs_sig$pval[glm_coeffs_sig$pval==0] = 2.2e-16
# glm_coeffs_sig$logPval = log10(glm_coeffs_sig$pval)
# glm_coeffs_sig$pval_direction = glm_coeffs_sig$logPval * ifelse(glm_coeffs_sig$coeff>=0,yes = 1, no = -1)
# ggplot(glm_coeffs_sig, aes(x = chromosome, y = feature, fill = pval_direction)) +
#   geom_raster() +
#   # scale_fill_gradient(low="grey90", high="red",na.value="grey") +
#   scale_fill_gradient2(
#     low = "blue",mid = "gray90",high = "red",
#     breaks = c(min(glm_coeffs_sig$pval_direction, na.rm = T), 0,
#                max(glm_coeffs_sig$pval_direction, na.rm = T)),
#     labels = c(round(min(glm_coeffs_sig$pval_direction, na.rm = T), digits = 1), 0,
#                round(max(glm_coeffs_sig$pval_direction, na.rm = T), digits = 1))) +
#   labs(y = "Predictor", x = "Chromosome") +
#   labs(fill = "log10 p-value")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2PWGS)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
#         axis.title=element_text(size=14,face="bold")) +
#   facet_wrap(c("tissue"), nrow = 1) +
#   geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")),
#             nudge_y = -0.4, col = "gray7", size=2)
# ggsave(paste0("fig/modelEvaluationWGS/GLMpvalsSig_perChr", plotEnding, ".png"),
#        height=8, width=13)
# ######
# 
# 
# # RF Predictor importances tissue comparison #####
# print("RF Predictor importances tissue comparison")
load("data/Modeling/WholeGenomeData/RF/rf_imps.RData")
rf_finalGini = rf_imps[rf_imps$chromosome == "final",]
ggplot(rf_finalGini, aes(x = tissue, y = feature, fill = gini_scaled)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  facet_grid(rows = vars(group), space = "free_y", scales = "free_y", switch = "y") +
  # labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Gini\nImportance")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2PWGS)+
  theme(axis.text.y=element_text(size=4.5),
        axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=4,face="bold"),
        strip.text.y.left = element_text(angle = 0, size=5.5),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines")) 
ggsave(paste0("fig/modelEvaluationWGS/RFginiScaled_Tissues", plotEnding, ".png"),
       height=9, width=6)
ggsave(paste0("fig/modelEvaluationWGS/RFginiScaled_Tissues", plotEnding, ".pdf"),
       height=9, width=6)

# #####
# 
# # GLM coeff pre-selection tissue comparison #####
# # not done in WGS
# #####
# 
# # GLM coeff post-selection tissue comparison #####
# print(" GLM coeff post-selection tissue comparison")
# glm_finalCoeffs_sig = glm_coeffs_sig[glm_coeffs_sig$chromosome == "final",]
# # coefficients
# ggplot(glm_finalCoeffs_sig, aes(x = tissue, y = feature, fill = coeff)) +
#   geom_raster() +
#   scale_fill_gradient2(
#     low = "blue",mid = "gray90",high = "red",na.value="grey",
#     breaks = c(min(glm_finalCoeffs_sig$coeff, na.rm = T), 0,
#                max(glm_finalCoeffs_sig$coeff, na.rm = T)),
#     labels = as.character(c(round(min(glm_finalCoeffs_sig$coeff, na.rm = T), digits = 1), 0,
#                             round(max(glm_finalCoeffs_sig$coeff, na.rm = T), digits = 1)))) +
#   labs(y = "Predictor", x = "Chromosome") +
#   labs(fill = "Coefficient")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2PWGS)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
#         axis.title=element_text(size=14,face="bold"))
# ggsave(paste0("fig/modelEvaluationWGS/GLMcoeffsSig_Tissues", plotEnding, ".png"),
#        height=8, width=6)
# # p-values
# ggplot(glm_finalCoeffs_sig, aes(x = tissue, y = feature, fill = pval_direction)) +
#   geom_raster() +
#   scale_fill_gradient2(
#     low = "blue",mid = "gray90",high = "red",
#     breaks = c(min(glm_finalCoeffs_sig$pval_direction, na.rm = T), 0,
#                max(glm_finalCoeffs_sig$pval_direction, na.rm = T)),
#     labels = c(round(min(glm_finalCoeffs_sig$pval_direction, na.rm = T), digits = 1), 0,
#                round(max(glm_finalCoeffs_sig$pval_direction, na.rm = T), digits = 1))) +
#   labs(y = "Predictor", x = "Chromosome") +
#   labs(fill = "log10 p-value")+
#   scale_x_discrete(labels=t2T) +
#   scale_y_discrete(labels = p2PWGS)+
#   theme(axis.text.y=element_text(size=5.5),
#         axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
#         axis.title=element_text(size=14,face="bold")) +
#   geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")),
#             nudge_y = -0.4, col = "gray7", size=2)
# ggsave(paste0("fig/modelEvaluationWGS/GLMpvalsSig_Tissues", plotEnding, ".png"),
#        height=8, width=6)
# #####



# direct comparison of importances #####
# print("direct comparison of importances")
# # rf_finalGini
# # glm_finalCoeffs_sig
# combinedCoeffs = data.frame("RF gini importance" = rf_finalGini$gini_scaled,
#                             "GLM coefficient" = glm_finalCoeffs_sig$coeff_scaled,
#                             "GLM log p-value" = glm_finalCoeffs_sig$logPval,
#                             tissue = rf_finalGini$tissue)
# pdf(paste0("fig/modelEvaluationWGS/compareMethodCoefficients",plotEnding,".pdf"))
# dumpVar = sapply(t2T, function(tissue){
#   plot(combinedCoeffs[combinedCoeffs$tissue == tissue,1:3], main = tissue)
# })
# dev.off()
#####


# TODO comparison WGS with WEX #####
#####