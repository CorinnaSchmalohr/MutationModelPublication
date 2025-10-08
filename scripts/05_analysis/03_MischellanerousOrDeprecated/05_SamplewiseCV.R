library(ROCR)
tissues = c("brain","breast", "colon","esophagus", 
            "kidney", "liver", "lung","ovary", 
            "prostate", "skin") #
t2T = setNames(c("Brain","Breast", "Colon","Esophagus", 
                 "Kidney", "Liver", "Lung","Ovary", 
                 "Prostate", "Skin"), tissues)
tissueCols = setNames(rainbow(length(tissues)), tissues)
darkCol = function(x,n){ # function to create a darkened version of a color
   sapply(x, function(y){
      cl = colorRampPalette(c(y, "black"))(n+1)[1:n]
   })
}
lightCol = function(x,alpha){
   x = col2rgb(x)
   rgb(red=x[1,], green=x[2,], blue=x[3,], alpha=alpha, maxColorValue=255)
}


# get the sample CV performances for RF and context #####
RFperfs= sapply(tissues, function(tissue){
   load(paste0("data/rdata/RFmodelPatientCV/", tissue,
               "_RFresults.RData"))
   preds = lapply(RFresults, function(x){x$predictions$pred})
   labels = lapply(RFresults, function(x){x$predictions$label})
   perf = prediction(preds, labels)
   roc = performance(perf, "tpr", "fpr")
   pr = performance(perf,"prec", "rec")
   auc = unlist(performance(perf,"auc")@y.values)
   return(list(roc = roc, pr = pr, auc = auc, samps = names(RFresults)))
}, simplify = F)

contextPerfs = sapply(tissues, function(tissue){
   print(tissue)
   perfs = sapply(RFperfs[[tissue]]$samps, function(samp){
      cat(samp,' ')
      preds = sapply(1:32, function(i){
         load(paste0("data/rdata/contextPatientCV/exomeMuts_part",
                     i,"_", tissue,"_",samp,"_results.RData")) 
         return(pred)
      }, simplify = F)
      preds = do.call(rbind,preds)
      perf =  prediction(preds$pred, preds$label)
      roc = performance(perf, "tpr", "fpr")
      pr = performance(perf,"prec", "rec")
      auc = performance(perf,"auc")@y.values[[1]]
      return(list(roc = roc, pr = pr, auc = auc))
   }, simplify=F); cat('\n')
   
   roc = perfs[[1]]$roc
   roc@x.values  = unname(lapply(perfs, function(x)x$roc@x.values[[1]]))
   roc@y.values = unname(lapply(perfs, function(x)x$roc@y.values[[1]]))
   roc@alpha.values = unname(lapply(perfs, function(x)x$roc@alpha.values[[1]]))
   
   pr = perfs[[1]]$pr
   pr@x.values  = unname(lapply(perfs, function(x)x$pr@x.values[[1]]))
   pr@y.values = unname(lapply(perfs, function(x)x$pr@y.values[[1]]))
   pr@alpha.values = unname(lapply(perfs, function(x)x$pr@alpha.values[[1]]))
   
   auc = sapply(perfs, function(x)x$auc)
   return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
save(contextPerfs, RFperfs, file = "data/rdata/patientCV_perfs.RData")
#####

# load #####
load("data/rdata/patientCV_perfs.RData")
# load the performances you get by using chromosome-wise CV
load("data/rdata/RFmodel/ROC_PR_RF_concat.RData") # ROC_PR_RF_concat
ROC_PR_context = sapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/procData/exomeMuts/combinedPredictions_",
               tissue, ".RData"))
   perf =  prediction(preds$contextPreds, preds$mutated)
   roc = performance(perf, "tpr", "fpr")
   pr = performance(perf,"prec", "rec")
   auc = performance(perf,"auc")@y.values[[1]]
   return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
load("data/rdata/RFmodel/RF_impsAndPvals.RData")

#####

# plot ROCs #####
png("fig/patienwiseCV_ROCs.png", width=1400, height = 1400, pointsize=30)
par(mfrow = c(5,4), mar = c(2,2,1,0), oma = c(1.5,1.5,0,0))
dumpVar = sapply(tissues, function(tissue){
   plot(RFperfs[[tissue]]$roc, col = "cornflowerblue", main = paste0(tissue, " RF"))
   plot(RFperfs[[tissue]]$roc, col = "blue", avg = "threshold", 
        add = T, lwd = 3)
   plot(ROC_PR_RF_concat[[tissue]]$roc, add = T, lwd = 3, col = "blue4",
        lty = 2)
   abline(0,1, lty = 2, col = "grey")
   plot(contextPerfs[[tissue]]$roc, col = "coral3", main = paste0(tissue, " context"))
   plot(contextPerfs[[tissue]]$roc, col = "red", add = T,
        avg = "threshold", lwd = 3)
   plot(ROC_PR_context[[tissue]]$roc, add = T, lwd = 3, col = "firebrick4", lty = 2)
   abline(0,1, lty = 2, col = "grey")
})
mtext("FPR", side = 1, outer = T)
mtext("TPR", side = 2, outer = T)
legend("bottomright", legend = c("patient-wise CV folds", 
                                 "average over patient-wise CV folds",
                                 "chromosome-wise CV"), 
       col = c( "grey", "grey", "black"), 
       lty = c(1,1,3), lwd = c(1,3,3))
dev.off()
#####


# plot PRs #####
png("fig/patienwiseCV_PRs.png", width=1400, height = 1400, pointsize=30)
par(mfrow = c(5,4), mar = c(2,2,1,0), oma = c(1.5,1.5,0,0))
dumpVar = sapply(tissues, function(tissue){
   plot(RFperfs[[tissue]]$pr, col = "cornflowerblue", main = paste0(tissue, " RF"))
   plot(RFperfs[[tissue]]$pr, col = "blue", avg = "threshold", 
        add = T, lwd = 3)
   plot(ROC_PR_RF_concat[[tissue]]$pr, add = T, lwd = 3, col = "blue4",
        lty = 2)
   plot(contextPerfs[[tissue]]$pr, col = "coral3", main = paste0(tissue, " context"))
   plot(contextPerfs[[tissue]]$pr, col = "red", add = T,
        avg = "threshold", lwd = 3)
   plot(ROC_PR_context[[tissue]]$pr, add = T, lwd = 3, col = "firebrick4", lty = 2)
})
mtext("Recall", side = 1, outer = T)
mtext("Precision", side = 2, outer = T)
legend("topright", legend = c("patient-wise CV folds", 
                                 "average over patient-wise CV folds",
                                 "chromosome-wise CV"), 
       col = c( "grey", "grey", "black"), 
       lty = c(1,1,3), lwd = c(1,3,3))
dev.off()
#####


# plot AUCs vs. data size ######
png("fig/patienwiseCV_AUCvsN.png", width=1000, height = 1400, pointsize=30)
par(mfrow = c(5,2), mar = c(3,3,1,0))
dumpVar = sapply(tissues, function(tissue){
   auc = RFperfs[[tissue]]$auc
   load(paste0("data/rdata/RFmodelPatientCV/", tissue,
               "_RFresults.RData"))
   ns = sapply(RFresults, function(x){nrow(x$predictions)})
   plot(ns, auc, main = tissue, mgp = c(2,0.8,0), xlab = "n positions", ylab = "AUC")
})
dev.off()
#####


# visualize predictor importances #####
load("data/rdata/RFmodel/RF_impsAndPvals.RData") #rf_imps,rf_gini

png("fig/patienwiseCV_importances.png", width=1000, height = 1400, pointsize=30)
par(mfrow = c(5,2), mar = c(2,2,1,0), oma = c(1.5,1.5,0,0))
dumpVar= sapply(tissues, function(tissue){
   print(tissue)
   load(paste0("data/rdata/RFmodelPatientCV/", tissue,
               "_RFresults.RData"))
   imps = sapply(RFresults, function(x){x$importance})
   impsM = rowMeans(imps)
   impsSD = apply(imps,1,sd)
   
   oldImps = rf_imps[[tissue]][rownames(imps),]
   oldM = rowMeans(oldImps)
   oldSD = apply(oldImps, 1,sd)
   
   plot(impsM, oldM,main = tissue)
   abline(0,1, col = "grey")
   segments(x0 = impsM-impsSD, x1 = impsM+impsSD, y0 = oldM, y1 = oldM)
   segments(x0 = impsM, x1 = impsM, y0 = oldM-oldSD, y1 = oldM+oldSD)
}, simplify = F)
title( xlab = "patient-wise CV permutation importance", 
       ylab ="chromosome-wise CV permutation importance", 
       outer = T, line = 0, cex.lab = 1.5)
dev.off()
#####

# combined plot only for RF #####
png("fig/patienwiseCV_Overview.png", width=1400, height = 1600, pointsize=30)
par(mfrow = c(10,4), mar = c(0.1,4,0,0), oma = c(3,1.5,0.2,0.2))
dumpVar = sapply(tissues, function(tissue){
   #roc
   plot(NA, xlim = c(0,1), ylim = c(0,1), las = 1, xaxt = "n", yaxt = "n",
        xlab = "", ylab = "")
   abline(0,1, lty = 2, col = "grey")
   axis(2,at = c(0,0.5,1), las = 1, mgp = c(2,0.7,0), padj=c(0.2,0.5,0.8))
   plot(RFperfs[[tissue]]$roc, col = lightCol(tissueCols[tissue], alpha = 70),
        add = T)
   plot(RFperfs[[tissue]]$roc, col = darkCol(tissueCols[tissue],5)[3], 
        avg = "threshold", add = T, lwd = 3)
   plot(ROC_PR_RF_concat[[tissue]]$roc, add = T, lwd = 3, 
        col = darkCol(tissueCols[tissue],4)[3],lty = 2)
   if(tissue == tail(tissues,1)){
      axis(1,at = c(0,0.5,1), las = 1, mgp = c(2,0.7,0))
      title(xlab = "FPR", xpd = NA, line = 1.5)
   }
   if(tissue == tissues[length(tissues)/2]){
      title(ylab = "TPR", xpd = NA, line = 2)
   }
   if(tissue == tissues[1]){
      legend("bottomright", legend = c("PWCV folds", 
                                       "aver. PWCV",
                                       "aver. CWCV"), 
             col = c( "grey", "dark grey", "black"), 
             lty = c(1,1,3), lwd = c(1,3,3), 
             y.intersp=0.8, seg.len=1.3, x.intersp=0.7)
   }
   mtext(t2T[tissue], side = 2, line = 4)
   # pr
   plot(NA, xlim = c(0,1), ylim = c(0,1), xaxt = "n", yaxt = "n",
        xlab = "", ylab = "")
   axis(2,at = c(0,0.5,1), las = 1, mgp = c(2,0.7,0), padj=c(0.2,0.5,0.8))
   plot(RFperfs[[tissue]]$pr, col = lightCol(tissueCols[tissue], alpha = 70),
        add = T)
   plot(RFperfs[[tissue]]$pr, col =  darkCol(tissueCols[tissue],5)[3], 
        avg = "threshold", add = T, lwd = 3)
   plot(ROC_PR_RF_concat[[tissue]]$pr, add = T, lwd = 3, 
        col = darkCol(tissueCols[tissue],4)[3],lty = 2)
   if(tissue == tail(tissues,1)){
      axis(1,at = c(0,0.5,1), las = 1, mgp = c(2,0.7,0))
      title(xlab = "Recall", xpd = NA, line = 1.5)
   }
   if(tissue == tissues[length(tissues)/2]){
      title(ylab = "Precision", xpd = NA, line = 2)
   }
   
   # auc vs. datasize
   auc = RFperfs[[tissue]]$auc
   load(paste0("data/rdata/RFmodelPatientCV/", tissue,
               "_RFresults.RData"))
   ns = sapply(RFresults, function(x){nrow(x$predictions)})
   plot( auc,  ns,mgp = c(2,0.8,0),  las = 1, xlab = "", ylab = "", xaxt = "n",
        xlim = range(sapply(RFperfs, function(x)range(x$auc))), log = "y")
   if(tissue == tail(tissues,1)){
      axis(1, las = 1, mgp = c(2,0.7,0))
      title(xlab = "AUC", xpd = NA, line = 1.5)
   }
   if(tissue == tissues[length(tissues)/2]){
      title(ylab = "Data size", xpd = NA, line = 3)
   }
   
   
   
   # imps
   load(paste0("data/rdata/RFmodelPatientCV/", tissue,
               "_RFresults.RData"))
   imps = sapply(RFresults, function(x){x$importance})
   impsM = rowMeans(imps)
   impsSD = apply(imps,1,sd)
   
   oldImps = rf_imps[[tissue]][rownames(imps),]
   oldM = rowMeans(oldImps)
   oldSD = apply(oldImps, 1,sd)
   
   plot(oldM, impsM, xlab = "", ylab = "", xaxt = "n", las = 1,
        xlim = range(sapply(rf_imps, range, na.rm = T)),
        ylim = range(sapply(rf_imps, range, na.rm = T)))
   abline(0,1, col = "grey")
   segments(y0 = impsM-impsSD, y1 = impsM+impsSD, x0 = oldM, x1 = oldM)
   segments(y0 = impsM, y1 = impsM, x0 = oldM-oldSD, x1 = oldM+oldSD)
   if(tissue == tail(tissues,1)){
      axis(1, las = 1, mgp = c(2,0.7,0))
      title(xlab = "CWCV importance", xpd = NA, line = 1.5)
   }
   if(tissue == tissues[length(tissues)/2]){
      title(ylab = "PWCV importance", xpd = NA, line = 3)
   }
})
dev.off()


#####
