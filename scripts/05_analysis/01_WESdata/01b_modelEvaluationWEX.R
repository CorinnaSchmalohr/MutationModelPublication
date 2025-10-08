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
dir.create("fig/modelEvaluation", showWarnings=F)
plotEnding = "_20230620"

# load data #####
load("data/processedData/dataInfos.RData")
load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
load("data/Modeling/exomeTrainData/GLM/predPerTissueGLM_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/predPerTissueLasso.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_perChr.RData")
load("data/Modeling/exomeTrainData/GLM/ROC_PR_glm_perChr_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_perChr.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
load("data/Modeling/exomeTrainData/GLM/ROC_PR_glm_concat_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/ROC_PR_lasso_concat.RData")
load("data/processedData/predictorOrder.RData")
load("data/Modeling/exomeTrainData/RF/RF_imps.RData")
load("data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals.RData")
load("data/Modeling/exomeTrainData/GLM/GLM_impsAndPvals_sig.RData")
load("data/Modeling/exomeTrainData/Lasso/Lasso_impAndStab.RData")
load("data/MutTables/exomeTrainData/nPerTissue.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_train.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_glm_sig_train.RData")
load("data/Modeling/exomeTrainData/RF/ROC_PR_lasso_sig_train.RData")
#####



# Corrplots ######
dumpVar = sapply(tissues, function(tissue){
  png(paste0("fig/modelEvaluation/corrplot_", tissue, plotEnding, ".png"),
      width=1550, height=1300, pointsize = 25)
  cors = dataInfos[[tissue]]$cors
  rownames(cors) = p2P[rownames(cors)]
  colnames(cors) = p2P[colnames(cors)]
  corrplot(cors, type = 'lower',  tl.col = 'black',      
           cl.ratio = 0.2, tl.cex=0.4, tl.srt = 45, col = COL2('PuOr', 10))
  dev.off()
})
#####

# proportion of muts removed by removing positions mutated more than once #####
load("data/MutTables/exomeTrainData/muts.RData")
tissue2Cancer = list("lung" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma",
                     "esophagus" = "Esophageal carcinoma",
                     "liver" = "Liver hepatocellular carcinoma",
                     "brain" = "Brain Lower Grade Glioma")
duplRemoved = sapply(tissues, function(tissue){
  print(tissue)
  # subset to tissue
  submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
  
  # exclude positions that were mutated more than once, 
  # to exclude possible selection
  positions = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
  multiple = duplicated(positions) | duplicated(positions, fromLast=T)
  table(multiple)
})
propRemoved = duplRemoved[2,]/colSums(duplRemoved)

pdfAndPng(file = "fig/duplicatedMutationsRemovedWEX", width = 10, height = 3, 
          expr = ({
            par(mar = c(3,5,1,1))
            temp = barplot(duplRemoved, las = 1, names.arg = t2T[tissues], 
                           legend.text = c("retained", "removed"), args.legend = list(x="topleft", inset = 0.01),
                           ylab = "n Mutations", mgp = c(4,1,0))
            text(x = temp, y = colSums(duplRemoved)+5000, 
                 labels = paste0(round(propRemoved*100, digits = 2), "%"), adj = c(0.5,0), xpd = NA)
          }))
#####

# Mutation overview #####
# prepare for upper part: mutation counts per sample, ordered
sortByMedian = nPerTissue[order(sapply(nPerTissue,median))]
# prepare for lower part: counts of mutation types
mutClass = list(c("T>A", "A>T"),
                c("T>C", "A>G"),
                c("T>G", "A>C"),
                c("C>A", "G>T"),
                c("C>G", "G>C"),
                c("C>T", "G>A"))
names(mutClass) = sapply(mutClass, paste, collapse = "/")
mutTypes = sapply(names(sortByMedian), function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped.RData"))
  muts = data$muts[data$muts$mutated == 1,]
  temp = paste(muts$ref, muts$alt, sep=">")
  sapply(mutClass, function(x){
    sum(temp %in% x)
  })
})
bases = c("A", "C", "G", "T")
trimClass = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
                                c(rep("C", 48), rep("T", 48)), ">",
                                c(rep(c("A", "G", "T"), each =16), 
                                  rep(c("A", "C", "G"), each =16)), "]",
                                rep(bases, 24))),
                  paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
                                c(rep("G", 48), rep("A", 48)), ">",
                                c(rep(c("T", "C", "A"), each =16), 
                                  rep(c("T", "G", "C"), each =16)), "]",
                                rep(rev(bases), 24))))
rownames(trimClass) = paste(trimClass[,1], trimClass[,2], sep = " and ")
trimTypes = sapply(names(sortByMedian), function(tissue){
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped.RData"))
  muts = data$muts[data$muts$mutated == 1,]
  temp = paste0(substr(muts$context,2,2),"[",muts$ref,">", muts$alt,"]",
                substr(muts$context,4,4))
  apply(trimClass,1, function(x){
    sum(temp %in% x)
  })
})

# plot
pdfAndPng(file = paste0("fig/MutationOverview", plotEnding),
          width = 29, height = 20, 
          pngArgs = list(res=200, pointsize=20),
          pdfArgs = list(pointsize=35), expr = expr({
            layout(rbind(c(1,1,1,1,1,1,1,1,1,1),
                         c(2,3,4,5,6,7,8,9,10,11),
                         c(2,3,4,5,6,7,8,9,10,11),
                         c(2,3,4,5,6,7,8,9,10,11),
                         c(2,3,4,5,6,7,8,9,10,11),
                         c(2,3,4,5,6,7,8,9,10,11)))
            # par(mar = c(0.1,0.15,0,0), oma = c(3.3,5,2,0))
            par(mar = c(0.1,1,0,0), oma = c(3.3,5,2,0))
            # upper part: mutation counts per sample, ordered
            sortByMedian = nPerTissue[order(sapply(nPerTissue,median))]
            plot(NA, xlim = c(1,20), ylim = c(1,max(unlist(nPerTissue))), log = "y",bty="l", 
                 ylab = "Mutations per tumor", xaxt = "n", las = 1, xlab = "", xpd = NA,
                 mgp = c(5,1,0), cex.lab = 1.1,cex.axis = 1.1)
            dumpVar = sapply(1:length(sortByMedian), function(i){
              tissue= names(sortByMedian)[i]
              x = sortByMedian[[i]]
              x = sort(x)
              pos = seq(0,1,length.out = length(x))
              points(i*2-1+pos,x, col = tissueCols[tissue], pch = 1, cex=0.7)
              segments(x0=i*2-1,x1=i*2, y0 = median(x), lwd = 2, lty = 2)
              # return(cbind(i*2-1+pos,x))
            })
            axis(3,at =1:length(sortByMedian)*2-0.5, labels= t2T[names(sortByMedian)], tick=F,
                 mgp = c(2,0.5,0), cex.axis  = 1.1)
            # lower part: frequencies of mutation types
            par(mar = c(0,1,0,0))
            temp =apply(trimTypes, 2,function(x){x/sum(x)})
            dumpVar = sapply(1:ncol(temp), function(i){
              barplot(rev(temp[,i]), horiz=T, las = 1,  space = 0, yaxs = "i",mgp = c(2,0.5,0),
                      names.arg="",cex.names=0.5, cex.axis = 1.1,
                      col = rep(rev(basesCol), each=16), xlim = c(0,0.16), xaxp = c(0, 0.1, 1))
              if(i==1){
                classes = trimClass[rownames(temp),1]
                names = rev(paste0(substr(classes,1,1), substr(classes,3,3), substr(classes,7,7)))
                names2 = rev(paste0(" ", substr(classes,3,3), " "))
                text(x = 0, y = 0.99:96-0.5,  labels=names,family = "mono",
                     xpd = NA, adj=1.1,  cex = 0.6)
                text(x = 0, y = 0.99:96-0.5,  labels=names2,family = "mono",
                     xpd = NA, adj=1.1, col = rep(rev(basesCol), each=16), font = 2, cex = 0.6)
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
          }))
#####

# Feature effect directions ######
preds = names(p2P)
for(tissue in tissues){
  print(tissue)
  load(paste0("data/MutTables/exomeTrainData/", 
              tissue, "_Muts_mapped_processed.RData"))
  pdf(paste0("fig/predictor_effectDirection_", tissue, plotEnding, ".pdf"),
      width = 20,height = 15,pointsize=20)
  par(mfrow = c(5,8), mar = c(2.5,3,3,0.5))
  sapply(preds, function(x){
    title = strwrap(p2P[x], width = 18)
    if(! x %in% colnames(dat)){
      # plot(0,type='n',axes=FALSE, main = title, xlab = "", ylab = "", cex.main=0.8)
    }else{
      nLevels = length(unique(dat[,x]))
      if(nLevels == 2){
        pd = cbind(table(dat[,x][dat$mutated == 0])/sum(dat$mutated == 0), NA)
        barplot(pd, names.arg=c("TN", "TP"),
                las = 1, main = title, cex.main=0.8,
                col = c(rgb(238,44,44,100, maxColorValue = 255),
                        rgb(238,44,44, maxColorValue = 255)), legend.text =F)
        pd = cbind(NA, table(dat[,x][dat$mutated == 1])/sum(dat$mutated == 1))
        barplot(pd, col = c(rgb(0,0,139,100, maxColorValue = 255),
                            rgb(0,0,139, maxColorValue = 255)),
                add= T, yaxt = "n", width=c(1,1))
      } else if (nLevels <=5){
        barplot(cbind(table(dat[,x][dat$mutated == 0]),
                      table(dat[,x][dat$mutated == 1]), NA), 
                names.arg=c("TN", "TP", ""),
                las = 1, main = title, cex.main=0.8,
                col = rainbow(nLevels), legend.text=T,
                args.legend = list(bty = "n"))
      }
      else if(is.numeric(dat[,x])){
        temp = density(dat[,x][dat$mutated == 1])
        temp2 = density(dat[,x][dat$mutated == 0])
        plot(temp, col = rgb(238,44,44, maxColorValue = 255), lwd = 2.5,
             main = title, xlab = "", las = 1,cex.main=0.8,
             ylim = c(0,max(c(temp$y, temp2$y))))
        polygon(temp, col = rgb(238,44,44,alpha = 75, maxColorValue = 255), border = NA)
        lines(temp2, col = rgb(0,0,139, maxColorValue = 255), lwd = 2.5, lty = 2)
        polygon(temp2, col = rgb(0,0,139,alpha = 75, maxColorValue = 255), border = NA)
      }  else
        stop("something went wrong")
    }
  })
  dev.off()
}
#####


# Comparison of RF, GLM and LASSO detailed plot #####
png(paste0("fig/modelEvaluation/compareMethodperformance", plotEnding,".png"), 
    height=1600, width=1300, pointsize=30)#, res = 200)
par(mfrow = c(length(tissues),3), mar = c(0.5,6,0.1,0.1), 
    mgp = c(2.5,1,0), oma = c(3,6,0.05,0.05))
lwd = 2
plotDump = sapply(tissues, function(tissue){
  print(tissue)
  # roc
  # plot(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]],
  #      ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]],lwd = 2, 
  #      col = methodCols["RF"], las = 1,xlab = "", 
  #      ylab = "", type = "l", xaxt = "n",yaxt = "n")
  # axis(2, mgp = c(2,0.7,0), las = 1,
  #      at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  # axis(2, mgp = c(2,0.7,0), las = 1,
  #      at = 1, padj = 0.7, cex.axis = 1.4)
  # abline(0,1, col = "grey", lty = 2, lwd = lwd)
  # plot(ROC_PR_glm_concat_sig[[tissue]]$roc, add = T,
  #      lwd = lwd, col = methodCols["GLM"], lty = 2)
  # plot(ROC_PR_lasso_concat[[tissue]]$roc, add = T,
  #      lwd = lwd, col = methodCols["SL"], lty = "4414")
  # add roc axis labels
  # if(tissue == tail(tissues,1)){   
  #   axis(1, mgp = c(2,0.7,0), las = 1,
  #        at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis=1.4)
  #   axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis=1.4)
  #   title(xlab = "FPR", line = 2, xpd = NA, cex.lab = 1.4)
  # }
  # if(tissue == tissues[ceiling(length(tissues)/2)]){
  #   title(ylab = "TPR", xpd = NA, mgp = c(2.5,1,0), cex.lab = 1.4)
  # }
  # auroc
  AUCs = data.frame(RF = sapply(ROC_PR_RF_perChr[[tissue]],function(x){x$auc}),
                    GLM = sapply(ROC_PR_glm_perChr_sig[[tissue]],function(x){x$auc}),
                    SL = sapply(ROC_PR_lasso_perChr[[tissue]],function(x){x$auc}))
  # boxplot(AUCs, las = 1, outwex = 1.5,
  #         border = methodCols, xaxt = "n", yaxt = "n")
  sinaplot(AUCs, adjust =0.1, maxwidth = 0.5, col = methodCols, xaxt = "n", yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,cex.axis = 1.4,
       at = par("yaxp")[1:2], lwd = 0, lwd.ticks=1)
  AUCs = c(RF = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
           GLM = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]],
           SL = ROC_PR_lasso_concat[[tissue]]$auc@y.values[[1]])
  points(AUCs, col = methodCols, pch = 19, cex = 2)
  mtext(side=2, text=t2T[tissue], line = 5, las = 1)
  # add auroc axis labels
  if(tissue == tail(tissues,1)){   
    # mtext(text=names(methodCols), side=1,  at=1:3, cex.lab = 1.4)
    axis(1,at = 1:3, labels=names(methodCols), mgp = c(2,0.7,0),
         cex.axis=1.4)
  }
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    # title(ylab = "AUC", xpd = NA, mgp = c(3,1,0), cex.lab = 1.4)
    mtext(text = "AUC", side = 2, adj = -0.5,cex.lab = 1.4, line = 3)
  }
  #prc
  plot(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
       ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], lwd = 3, 
       las = 1,xlab = "", xaxt = "n", type = "l", col= methodCols["RF"],
       ylab = "",  ylim = c(0,1), mgp =  c(2,0.7,0),  xaxt = "n",yaxt = "n")
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  plot(ROC_PR_glm_concat_sig[[tissue]]$pr, lwd = 3, 
       col = methodCols["GLM"], add = T, lty = 2)
  plot(ROC_PR_lasso_concat[[tissue]]$pr, lwd = 3, 
       col = methodCols["SL"], add = T, lty = "4414")
  rect(xleft = -0.005, xright=0.02, ybottom=0.5, ytop=1, lwd = 1.5)
  # add pr axis labels
  if(tissue == tail(tissues,1)){   
    axis(1, mgp = c(2,0.7,0), las = 1,
         at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
    axis(1, mgp = c(2,0.7,0), las = 1,at = 1, cex.axis = 1.4)
    title(xlab = "Recall", line = 2, xpd = NA, cex.lab = 1.4)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    mtext(text = "Precision", side = 2, adj = -5,cex.lab = 1.4, line = 3)
  }

  # create inset with zoom
  lim = par("plt")
  xlims = lim[2]-lim[1]
  ylims = lim[4]-lim[3]
  inset = 0.05
  smallPlot(expr={
    plot(x = ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]],
         y = ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]],
         xlab = "", ylab = "", xaxt = "n", yaxt = "n",
         xlim = c(0,0.02),ylim = c(0.5,1), 
         lwd = 3, type = "l", col = methodCols["RF"])
    lines(x = ROC_PR_glm_concat_sig[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_glm_concat_sig[[tissue]]$pr@y.values[[1]],
          col = methodCols["GLM"], lty = 2, lwd = 3)
    lines(x = ROC_PR_lasso_concat[[tissue]]$pr@x.values[[1]],
          y = ROC_PR_lasso_concat[[tissue]]$pr@y.values[[1]],
          col = methodCols["SL"], lty = 2, lwd = 3)
    box(lwd = 2)},  
    x1 = lim[1]+0.3*xlims, x2 = lim[2]-0.05*xlims,
    y1 = lim[3]+0.05*ylims, y2 = lim[4]-0.3*ylims, xpd = F,
    mar = c(0,0,0,0), border = "transparent")
  
  
  # violin of predictions
  rfpreds = do.call(rbind,predPerTissueRF[[tissue]])
  glmpreds = do.call(rbind,predPerTissueGLMsig[[tissue]])
  lassopreds =  do.call(rbind,predPerTissueLasso[[tissue]])
  rfpredssplit = split(rfpreds$pred, rfpreds$label)
  glmpredssplit = split(glmpreds$pred, glmpreds$label)
  lassopredssplit = split(lassopreds$pred, lassopreds$label)
  
  plotDat = list("TN.RF" = rfpredssplit$`0`,
                 "TP.RF" = rfpredssplit$`1`,
                 "TN.GLM" = glmpredssplit$`0`,
                 "TP.GLM" = glmpredssplit$`1`,
                 "TN.SL" = lassopredssplit$`0`,
                 "TP.SL" = lassopredssplit$`1`)
  tempCols = RColorBrewer::brewer.pal(3,"Set2")
  sinaplot(plotDat, col = c(tempCols,methodCols)[c(1,4,2,5,3,6)], 
           xaxt = "n",yaxt = "n",las = 1,  ylab = "", cex = 0.8, ylim= c(0,1))
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = c(0,0.5), lwd = 0, lwd.ticks=1, cex.axis = 1.4)
  axis(2, mgp = c(2,0.7,0), las = 1,
       at = 1, padj = 0.7, cex.axis = 1.4)
  abline(h=0.5, col = "grey", lty = 2, lwd = 2)
  abline(v=c(2.5,4.5))
  boxplot(plotDat, col = c("grey80", "grey40"),
          add = T, boxwex = 0.2, outline = F, #col = rgb(0,0,0,alpha = 0),
          ann = F, yaxt = "n", xaxt = "n", mgp = c(2,0.7,0))
  if(tissue == tail(tissues,1)){   
    axis(1,at = 1:6, labels=rep(c("0","1"),3), mgp = c(2,0.7,0),
         cex.axis=1.2)
    title(xlab = "True labels", line = 2, xpd = NA, 
          mgp = c(2,0.7,0), cex.lab = 1.4)
    # axis(1,at = c(1.5,3.5,5.5), labels=names(methodCols), mgp = c(2,0.7,0),
    #      cex.axis=1.2, tick = F)
    text(x= c(1.5,3.5,5.5), y=0.06,labels=names(methodCols), cex = 1.2, font = 2)
  }  
  if(tissue == tissues[ceiling(length(tissues)/2)]){
    mtext(text = "Prediction", side = 2, adj = -15, cex.lab = 1.4, line = 3)
  }
})
dev.off()   
#####


# Comparison of AUC between RF, GLM and LASSO #####
AUCcollection = sapply(tissues, function(tissue){
  c(RF = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]],
    GLM = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]],
    SL = ROC_PR_lasso_concat[[tissue]]$auc@y.values[[1]])
})
png(paste0("fig/modelEvaluation/AUCbetweenMethods_allTissues",
           plotEnding, ".png"), width = 1200, height = 500, pointsize = 20)
par(mar = c(3,4,1,4))
temp = barplot(AUCcollection-0.5, beside = T, las = 1, ylab = "AUC",
        col = methodCols, axisnames = F, yaxt = "n")
axis(2, at = axTicks(2), labels = axTicks(2)+0.5, las = 1)
title(ylab = "AUC")
legend(x = 41, y=mean(AUCcollection)-0.5, legend = names(methodCols), 
       fill = methodCols, xpd = NA, bty = "n")
mtext(t2T[tissues], side = 1, at = temp[2,])
dev.off()
#####

# Train vs. test performance of RF, GLM and LASSO #####
RFaucs = sapply(tissues, function(tissue){
  c(train = ROC_PR_RF_train[[tissue]]$auc@y.values[[1]], 
    test = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]])
})
GLMaucs = sapply(tissues, function(tissue){
  c(train = ROC_PR_glm_sig_train[[tissue]]$auc@y.values[[1]], 
    test = ROC_PR_glm_concat_sig[[tissue]]$auc@y.values[[1]])
})
LASSOaucs = sapply(tissues, function(tissue){
  c(train = ROC_PR_lasso_sig_train[[tissue]]$auc@y.values[[1]], 
    test = ROC_PR_lasso_concat[[tissue]]$auc@y.values[[1]])
})
png(paste0("fig/modelEvaluation/TrainVsTestAUC",
           plotEnding, ".png"), width = 1200, height = 500, pointsize = 25)
par(mfrow = c(3,1), mar = c(1,4,0.5,0.5), oma = c(2,4,2,0))
barplot(RFaucs, beside = T, las = 1,
        ylab = "AUC", xaxt = "n",
        col = rep(tissueCols, each = 2), 
        density = rep(c(NA, 20),length(tissues)),
        legend.text = c("Training set", "Test set"), 
        args.legend = list(x="topright", density=c(NA,40), ncol = 2,
                           inset=c(0.05,-0.43), xpd = NA, fill = "black"))
mtext(text = "RF", side = 2, las = 1,  line = 5,) 
abline(h=0)
barplot(GLMaucs, beside = T, las = 1,
         ylab = "AUC", xaxt = "n",
        col = rep(tissueCols, each = 2), 
        density = rep(c(NA, 20),length(tissues)))
mtext(text = "GLM", side = 2, las = 1,  line = 5) 
abline(h=0)
barplot(LASSOaucs, beside = T, las = 1,
         ylab = "AUC", names.arg = t2T[tissues],
        col = rep(tissueCols, each = 2), 
        density = rep(c(NA, 20),length(tissues)))
mtext(text = "SL", side = 2, las = 1,  line = 5) 
abline(h=0)
dev.off()
#####

# ROC over chrs for all tissues #####
png(paste0("fig/modelEvaluation/RF_ROCsOverChrsAllTissues",plotEnding,".png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.7,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  abline(0,1, col = "grey", lty = 2)
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_RF_perChr[[tissue]][[cr]]$roc, 
         col = rgb(0,0,0,0.5), add = T)
  })
  lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
        ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
        col = tissueCols[tissue], lwd = 2.5)
  text(x=0.05,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(ylab = "TPR", mgp = c(2,0.7,0), xpd = NA, outer=T)
title(xlab = "FPR", mgp = c(1.8,0.7,0), xpd = NA, outer=T)
legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
dev.off()
#####


# PR over chrs for all tissues #####
png(paste0("fig/modelEvaluation/RF_PRsOverChrsAllTissues", plotEnding,".png"),
    width=1200, height=400, res=150)
par(oma = c(3,3,0.2,0.2), mar = c(0.5,0.5,0,00), mfrow = c(2,5))
plotDump = sapply(tissues, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       mgp = c(2,0.8,0), las = 1,xlab = "", ylab = "",
       xaxt = "n", yaxt = "n")
  dump = sapply(names(chrCols), function(cr){
    plot(ROC_PR_RF_perChr[[tissue]][[cr]]$pr, 
         col = rgb(0,0,0,0.5), add = T)
  })
  plot(ROC_PR_RF_concat[[tissue]]$pr, 
       col = tissueCols[tissue], add = T, lwd = 2.5, main = "")
  text(x=0.5,y=0.9,  labels=t2T[tissue],font= 2, cex= 1.1, adj= c(0.5,0))
  if(tissue %in% c("brain", "liver")){
    axis(2, las = 1, mgp = c(2,0.7,0))
  }
  if(tissue %in% tissues[6:10]){
    axis(1, mgp = c(2,0.7,0))
  }
})
title(xlab = "Recall", mgp = c(1.8,0.7,0), xpd = NA, outer = T)
title(ylab = "Precision", mgp = c(2,0.7,0), xpd = NA, outer = T)
legend("bottomright", lty = 1, lwd = c(1,2.5), bty = "n",
       col = c(rgb(0,0,0,0.5), "black"), legend=c("CWCV", "Total"))
dev.off()
#####

# AUC versus chromosome features #####
load("data/processedData/codExons_noUTR_filtered.RData")
load("data/processedData/chrLengths.RData")
codingLength = sapply(split(codExons_filtered$width, 
                            codExons_filtered$chr), sum)
pdfAndPng(paste0("fig/modelEvaluation/RF_chromAUCsFeatures", plotEnding),
          width=12, height=4, pngArgs = list(pointsize=17), pdfArgs = list(pointsize=14), 
          expr = expression({
            par(mfrow = c(3,10), mar = c(2,2,1,0.5), oma = c(1,2.5,1,0.5))
            # AUC versus training data size
            dumpVar = sapply(tissues, function(tissue){
              AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
              nMuts = dataInfos[[tissue]]$nMutsPerChr/1000
              plot(nMuts, AUCs, las = 1, mgp = c(2,0.7,0),
                   ylab = "", xlab = "", main = t2T[tissue])
              mod = lm(AUCs~nMuts)
              abline(mod, col = "grey", lty = 2, lwd = 2)
              if(tissue == tissues[1]){
                mtext(text = "AUC", side = 2, line = 3)
              }
              if(tissue == tissues[length(tissues)/2]){
                mtext(text = "Chromsome data size (1000 positions)", 
                      side = 1, line = 1.8, cex = 0.8)
              }
            })
            # AUC versus chromosome coding region length
            xVar = codingLength
            dumpVar = sapply(tissues, function(tissue){
              AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
              plot(xVar/1000000, AUCs, las = 1, mgp = c(2,0.7,0),
                   ylab = "", xlab = "")
              mod = lm(AUCs~xVar)
              abline(mod, col = "grey", lty = 2, lwd = 2)
              if(tissue == tissues[1]){
                mtext(text = "AUC", side = 2, line = 3)
              }
              if(tissue == tissues[length(tissues)/2]){
                mtext(text = "Coding region length (Mb)", side = 1, line = 1.8, cex = 0.8)
              }
            })
            # AUC versus gene density/percent coding
            xVar = codingLength/chrLengths
            dumpVar = sapply(tissues, function(tissue){
              AUCs = sapply(ROC_PR_RF_perChr[[tissue]], function(y){y$auc})
              plot(xVar, AUCs, las = 1,mgp = c(2,0.7,0),
                   ylab = "", xlab = "")
              mod = lm(AUCs~xVar)
              abline(mod, col = "grey", lty = 2, lwd = 2)
              if(tissue == tissues[1]){
                mtext(text = "AUC", side = 2, line = 3)
              }
              if(tissue == tissues[length(tissues)/2]){
                mtext(text = "Percent coding", side = 1, line = 1.8, cex = 0.8)
              }
            })
          }))

#####


# Overview for tissue comparison #####
pdfAndPng(file = paste0("fig/modelEvaluation/compareRFperformanceTissues", plotEnding),
          width=24, height=8, pngArgs = list(pointsize=45),pdfArgs = list(pointsize=30), 
          expr = expr({
            par(mfrow = c(1,3), mar = c(4,4,0.1,0.1), mgp = c(2.5,1,0))
            # roc
            plot(NA, xlim = c(0,1), ylim = c(0,1), 
                 xlab = "FPR",ylab = "TPR", las = 1, cex.lab = 1.2)
            abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
            plotDump = sapply(tissues, function(tissue){
              lines(ROC_PR_RF_concat[[tissue]]$roc@x.values[[1]], 
                    ROC_PR_RF_concat[[tissue]]$roc@y.values[[1]], 
                    col = tissueCols[tissue], lwd = 2)
            })
            
            # pr
            plot(NA, xlim = c(0,1), ylim = c(0,1), 
                 xlab = "Recall", ylab = "Precision", las = 1, cex.lab = 1.2)
            plotDump = sapply(tissues, function(tissue){
              lines(ROC_PR_RF_concat[[tissue]]$pr@x.values[[1]], 
                    ROC_PR_RF_concat[[tissue]]$pr@y.values[[1]], 
                    col = tissueCols[tissue], lwd = 2)
            })
            legend("bottomright", col=tissueCols, legend=t2T[names(tissueCols)],
                   lwd = 2, ncol = 2)
            ## AUCs vs nMuts ##
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
            
            plot(nMuts, AUCs, las = 1,
                 bg = tissueCols, pch = 21, cex = 1.5, cex.lab = 1.2,
                 ylim = c(min(colMeans(AUCsPerChr)-AUCsSD), 
                          max(colMeans(AUCsPerChr)+AUCsSD)), 
                 xlim = c(min(nMuts), max(nMuts)+50000), 
                 xlab = "n Positions", ylab = "AUC")
            abline(lm(colMeans(AUCsPerChr) ~ nMuts), col = "grey", lty = 2)
            arrows(y0 = colMeans(AUCsPerChr)-AUCsSD, x0 = nMuts,
                   y1 = colMeans(AUCsPerChr)+AUCsSD, x1 = nMuts, code = 3,
                   angle = 90, length = 0.1)
            points(nMuts, AUCs,
                   bg = tissueCols, pch = 21, cex = 1.5)
            # points(x=nMuts, y = AUCs, pch = 4, bg = tissueCols, cex = 1.5)
            text(x=nMuts, y = AUCs, labels=t2T[tissues], pos=4, cex = 1.2)
          }))
# png(paste0("fig/modelEvaluation/compareRFperformanceTissues", plotEnding, ".png"), 
#     width=2400, height=800, pointsize=35)

dev.off()
#####


# Mbp-wise performance #####
# only makes sense for WGS
#####

# Predictor importances over chromosomes RF ######
# rf gini importance
rf_imps = do.call(rbind,sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_", 
         "finalModel.RData")) # rf, importance
  temp = cbind(rf_gini[[tissue]],final = importance[rownames(rf_gini[[tissue]])])
  temp =as.data.table(temp, keep.rownames = "feature")
  temp$feature = factor(temp$feature, levels = temp$feature)
  temp = melt(temp, id.vars = "feature", 
              variable.name = "chromosome", value.name = "gini")
  temp$tissue = t2T[tissue]
  temp$gini_scaled = scale(temp$gini, center = F)
  return(temp)
}, simplify=F))
rf_imps$group = factor(p2G[rf_imps$feature], levels = rev(unique(p2G)))
ggplot(rf_imps, aes(x = chromosome, y = feature, fill = gini_scaled)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Scaled Gini\nImportance")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5),
        axis.text.x=element_blank(),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1)
ggsave(paste0("fig/modelEvaluation/RFgini_perChr", plotEnding, ".png"), 
       height=13, width=13)
ggsave(paste0("fig/modelEvaluation/RFgini_perChr", plotEnding, ".pdf"), 
       height=13, width=13)
######

# Predictor importances over chromosomes GLM pre-selection ######
# coefficients
# collect coefficients, add info if significant
glm_coeffs= do.call(rbind,sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/GLM/", tissue, ".RData"))
  # coefficients
  temp = cbind(glm_imps[[tissue]],final = logR$coefficients[names(p2P)])
  temp = as.data.table(temp, keep.rownames = "feature")
  temp$feature = factor(temp$feature, levels = temp$feature)
  temp = melt(temp, id.vars = "feature", 
              variable.name = "chromosome", value.name = "coeff")
  # and pvalues
  temp2 = cbind(glm_pvals[[tissue]],final = coef(summary(logR))[,4][-1][names(p2P)])
  temp2 = as.data.table(temp2, keep.rownames = "feature")
  temp2$feature = factor(temp2$feature, levels = temp2$feature)
  temp2 = melt(temp2, id.vars = "feature", 
              variable.name = "chromosome", value.name = "pval")
  # add pvalues to coefficient table
  temp$pval = temp2$pval
  temp$tissue = t2T[tissue]
  temp$coeff_scaled = scale(temp$coeff, center = F)
  return(temp)
},simplify=F))
ggplot(glm_coeffs, aes(x = chromosome, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_coeffs$coeff, na.rm = T), 0,
               max(glm_coeffs$coeff, na.rm = T)),
    labels = c(round(min(glm_coeffs$coeff, na.rm = T), digits = -1), 0,
               round(max(glm_coeffs$coeff, na.rm = T), digits = -1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1)
ggsave(paste0("fig/modelEvaluation/GLMcoeff_perChr", plotEnding, ".png"), 
       height=8, width=13)
# pvalues
glm_coeffs$pval[glm_coeffs$pval==0] = 2.2e-16
glm_coeffs$logPval = log10(glm_coeffs$pval)
glm_coeffs$pval_direction = glm_coeffs$logPval * ifelse(glm_coeffs$coeff>=0,yes = 1, no = -1)
ggplot(glm_coeffs, aes(x = chromosome, y = feature, fill = pval_direction)) + 
  geom_raster() +
  # scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_coeffs$pval_direction, na.rm = T), 0,
               max(glm_coeffs$pval_direction, na.rm = T)),
    labels = c(round(min(glm_coeffs$pval_direction, na.rm = T), digits = 1), 0,
               round(max(glm_coeffs$pval_direction, na.rm = T), digits = 1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "log10 p-value")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1) +
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/GLMpvals_perChr", plotEnding, ".png"), 
       height=8, width=13)
######

# Predictor importances over chromosomes GLM post-selection ######
# coefficients
# collect coefficients, add info if significant
glm_coeffs_sig= do.call(rbind,sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/GLM/", tissue, "_sig.RData"))
  # coefficients
  temp = cbind(glm_imps_sig[[tissue]],final = logR$coefficients[names(p2P)])
  temp = as.data.table(temp, keep.rownames = "feature")
  temp$feature = factor(temp$feature, levels = temp$feature)
  temp = melt(temp, id.vars = "feature", 
              variable.name = "chromosome", value.name = "coeff")
  # and pvalues
  temp2 = cbind(glm_pvals_sig[[tissue]],final = coef(summary(logR))[,4][-1][names(p2P)])
  temp2 = as.data.table(temp2, keep.rownames = "feature")
  temp2$feature = factor(temp2$feature, levels = temp2$feature)
  temp2 = melt(temp2, id.vars = "feature", 
               variable.name = "chromosome", value.name = "pval")
  # add pvalues to coefficient table
  temp$pval = temp2$pval
  temp$tissue = t2T[tissue]
  temp$coeff_scaled = scale(temp$coeff, center = F)
  return(temp)
},simplify=F))
ggplot(glm_coeffs_sig, aes(x = chromosome, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_coeffs_sig$coeff, na.rm = T), 0,
               max(glm_coeffs_sig$coeff, na.rm = T)),
    labels = c(round(min(glm_coeffs_sig$coeff, na.rm = T), digits = -1), 0,
               round(max(glm_coeffs_sig$coeff, na.rm = T), digits = -1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1)
ggsave(paste0("fig/modelEvaluation/GLMcoeffSig_perChr", plotEnding, ".png"), 
       height=8, width=13)
# pvalues
glm_coeffs_sig$pval[glm_coeffs_sig$pval==0] = 2.2e-16
glm_coeffs_sig$logPval = log10(glm_coeffs_sig$pval)
glm_coeffs_sig$pval_direction = glm_coeffs_sig$logPval * ifelse(glm_coeffs_sig$coeff>=0,yes = 1, no = -1)
ggplot(glm_coeffs_sig, aes(x = chromosome, y = feature, fill = pval_direction)) + 
  geom_raster() +
  # scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_coeffs_sig$pval_direction, na.rm = T), 0,
               max(glm_coeffs_sig$pval_direction, na.rm = T)),
    labels = c(round(min(glm_coeffs_sig$pval_direction, na.rm = T), digits = 1), 0,
               round(max(glm_coeffs_sig$pval_direction, na.rm = T), digits = 1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "log10 p-value")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1) +
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/GLMpvalsSig_perChr", plotEnding, ".png"), 
       height=8, width=13)
######

# Predictor importances over chromosomes LASSO ######
lasso_coeffs= do.call(rbind,sapply(tissues, function(tissue){
  load(paste0("data/Modeling/exomeTrainData/Lasso/", tissue, ".RData")) # sp, stab
  stability = sp$x[,stab$lpos][names(p2P)]
  load(paste0("data/Modeling/exomeTrainData/Lasso/", tissue, "_sig.RData")) #logR, sigFeatures
  importance = logR$coefficients[names(p2P)]
  
  # coefficients
  temp = cbind(lasso_imp[[tissue]],final = logR$coefficients[names(p2P)])
  temp = as.data.table(temp, keep.rownames = "feature")
  temp$feature = factor(temp$feature, levels = temp$feature)
  temp = melt(temp, id.vars = "feature", 
              variable.name = "chromosome", value.name = "coeff")
  # and stability
  temp2 = cbind(lasso_stability[[tissue]],final = sp$x[,stab$lpos][names(p2P)])
  temp2 = as.data.table(temp2, keep.rownames = "feature")
  temp2$feature = factor(temp2$feature, levels = temp2$feature)
  temp2 = melt(temp2, id.vars = "feature", 
               variable.name = "chromosome", value.name = "stab")
  # add pvalues to coefficient table
  temp$stab = temp2$stab
  temp$tissue = t2T[tissue]
  temp$coeff_scaled = scale(temp$coeff, center = F)
  return(temp)
}, simplify=F))
#  stability
ggplot(lasso_coeffs, aes(x = chromosome, y = feature, fill = stab)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Stability")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1) +
  geom_text(aes(label = ifelse(stab >= 0.6, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/LassoStability_perChr", plotEnding, ".png"), 
       height=8, width=13)
#  importance
ggplot(lasso_coeffs, aes(x = chromosome, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(lasso_coeffs$coeff, na.rm = T), 
               max(lasso_coeffs$coeff, na.rm = T)),
    labels = c(round(min(lasso_coeffs$coeff, na.rm = T), digits = -1), 
               round(max(lasso_coeffs$coeff, na.rm = T), digits = -1))) +
  # scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=5.5 , angle = 90,vjust = 0.5, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  facet_wrap(c("tissue"), nrow = 1) 
ggsave(paste0("fig/modelEvaluation/LassoCoeffs_perChr", plotEnding, ".png"), 
       height=8, width=13)
######

# RF Predictor importances tissue comparison #####
rf_finalGini = rf_imps[rf_imps$chromosome == "final",]
ggplot(rf_finalGini, aes(x = tissue, y = feature, fill = gini)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  facet_grid(rows = vars(group), space = "free_y", scales = "free_y", switch = "y") +
  # labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Gini\nImportance")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=4.5),
        axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=4,face="bold"),
        strip.text.y.left = element_text(angle = 0, size=5.5),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines")) 
ggsave(paste0("fig/modelEvaluation/RFgini_Tissues", plotEnding, ".png"), 
       height=8, width=6)
ggsave(paste0("fig/modelEvaluation/RFgini_Tissues", plotEnding, ".pdf"), 
       height=8, width=6)
ggplot(rf_finalGini, aes(x = tissue, y = feature, fill = gini_scaled)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  facet_grid(rows = vars(group), space = "free_y", scales = "free_y", switch = "y") +
  # labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Gini\nImportance")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=4.5),
        axis.text.x=element_text(size=8 , angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=4,face="bold"),
        strip.text.y.left = element_text(angle = 0, size=5.5),
        strip.placement = "outside",
        panel.spacing = unit(0.1, "lines")) 
ggsave(paste0("fig/modelEvaluation/RFginiScaled_Tissues", plotEnding, ".png"), 
       height=8, width=6)
ggsave(paste0("fig/modelEvaluation/RFginiScaled_Tissues", plotEnding, ".pdf"), 
       height=8, width=6)
#####

# GLM coeff pre-selection tissue comparison #####
glm_finalCoeffs = glm_coeffs[glm_coeffs$chromosome == "final",]
# coefficients
ggplot(glm_finalCoeffs, aes(x = tissue, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red",na.value="grey",
    breaks = c(min(glm_finalCoeffs$coeff, na.rm = T), 0,
               max(glm_finalCoeffs$coeff, na.rm = T)),
    labels = as.character(c(round(min(glm_finalCoeffs$coeff, na.rm = T), digits = 1), 0,
               round(max(glm_finalCoeffs$coeff, na.rm = T), digits = 1)))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold")) 
ggsave(paste0("fig/modelEvaluation/GLMcoeffs_Tissues", plotEnding, ".png"), 
       height=8, width=6)
# p-values
ggplot(glm_finalCoeffs, aes(x = tissue, y = feature, fill = pval_direction)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_finalCoeffs$pval_direction, na.rm = T), 0,
               max(glm_finalCoeffs$pval_direction, na.rm = T)),
    labels = c(round(min(glm_finalCoeffs$pval_direction, na.rm = T), digits = 1), 0,
               round(max(glm_finalCoeffs$pval_direction, na.rm = T), digits = 1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "log10 p-value")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/GLMpvals_Tissues", plotEnding, ".png"), 
       height=8, width=6)
#####

# GLM coeff post-selection tissue comparison #####
glm_finalCoeffs_sig = glm_coeffs_sig[glm_coeffs_sig$chromosome == "final",]
# coefficients
ggplot(glm_finalCoeffs_sig, aes(x = tissue, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red",na.value="grey",
    breaks = c(min(glm_finalCoeffs_sig$coeff, na.rm = T), 0,
               max(glm_finalCoeffs_sig$coeff, na.rm = T)),
    labels = as.character(c(round(min(glm_finalCoeffs_sig$coeff, na.rm = T), digits = 1), 0,
                            round(max(glm_finalCoeffs_sig$coeff, na.rm = T), digits = 1)))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold")) 
ggsave(paste0("fig/modelEvaluation/GLMcoeffsSig_Tissues", plotEnding, ".png"), 
       height=8, width=6)
# p-values
ggplot(glm_finalCoeffs_sig, aes(x = tissue, y = feature, fill = pval_direction)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(glm_finalCoeffs_sig$pval_direction, na.rm = T), 0,
               max(glm_finalCoeffs_sig$pval_direction, na.rm = T)),
    labels = c(round(min(glm_finalCoeffs_sig$pval_direction, na.rm = T), digits = 1), 0,
               round(max(glm_finalCoeffs_sig$pval_direction, na.rm = T), digits = 1))) +
  labs(y = "Predictor", x = "Chromosome") + 
  labs(fill = "log10 p-value")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold")) +
  geom_text(aes(label = ifelse(pval < 0.05, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/GLMpvalsSig_Tissues", plotEnding, ".png"), 
       height=8, width=6)
#####

# LASSO importances tissue comparison #####
lasso_finalCoeffs = lasso_coeffs[lasso_coeffs$chromosome == "final",]
ggplot(lasso_finalCoeffs, aes(x = tissue, y = feature, fill = stab)) + 
  geom_raster() +
  scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Stability")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold"))  +
  geom_text(aes(label = ifelse(stab >= 0.6, yes = "*", no = " ")), 
            nudge_y = -0.4, col = "gray7", size=2)
ggsave(paste0("fig/modelEvaluation/LassoStability_Tissues", plotEnding, ".png"), 
       height=8, width=6)
ggplot(lasso_finalCoeffs, aes(x = tissue, y = feature, fill = coeff)) + 
  geom_raster() +
  scale_fill_gradient2(
    low = "blue",mid = "gray90",high = "red", 
    breaks = c(min(lasso_finalCoeffs$coeff, na.rm = T), 
               max(lasso_finalCoeffs$coeff, na.rm = T)),
    labels = c(round(min(lasso_finalCoeffs$coeff, na.rm = T), digits = -1), 
               round(max(lasso_finalCoeffs$coeff, na.rm = T), digits = -1))) +
  # scale_fill_gradient(low="grey90", high="red",na.value="grey") +
  labs(y = "Predictor", x = "Tissue") + 
  labs(fill = "Coefficient")+
  scale_x_discrete(labels=t2T) +
  scale_y_discrete(labels = p2P)+
  theme(axis.text.y=element_text(size=5.5),
        axis.text.x=element_text(size=8, angle = 45,vjust = 1, hjust=1),
        axis.title=element_text(size=14,face="bold")) 
ggsave(paste0("fig/modelEvaluation/LassoCoeffs_Tissues", plotEnding, ".png"), 
       height=8, width=6)
#####


# direct comparison of importances #####
rf_finalGini
glm_finalCoeffs_sig
lasso_finalCoeffs
if(!all.equal(rf_finalGini$feature, glm_finalCoeffs$feature) |
   !all.equal(rf_finalGini$feature, lasso_finalCoeffs$feature) |
   !all.equal(rf_finalGini$tissue, lasso_finalCoeffs$tissue)){
  print("WARNING!")
}
combinedCoeffs = data.frame("RF gini importance" = rf_finalGini$gini_scaled, 
                        "GLM coefficient" = glm_finalCoeffs$coeff_scaled, 
                       "GLM log p-value" = glm_finalCoeffs$logPval,
                       "SL coefficient" = lasso_finalCoeffs$coeff_scaled, 
                       "SL stability" = lasso_finalCoeffs$stab, 
                       tissue = rf_finalGini$tissue)
pdf(paste0("fig/modelEvaluation/compareMethodCoefficients",plotEnding,".pdf"))
dumpVar = sapply(t2T, function(tissue){
  plot(combinedCoeffs[combinedCoeffs$tissue == tissue,1:5], main = tissue)
})
dev.off()
#####