# preparation #####
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
library(ROCR)
library("RColorBrewer")
library(viridis)
#####


# collect predictions for RF and context ######
print("collecting predictions")
dumpVar = sapply(tissues, function(tissue){
  print(tissue)
  # load information which positions were mutated
  load(paste0("data/procData/traindata/traindata_", tissue, ".RData"))
  tMuts = data$muts
  #toExclude = data$pred$ConsensusExcludable == 1 | data$pred$repeatMasker == 1 | data$pred$tandemRepeatFinder == 1 | 
  #  data$pred$mappability_100mer < 1 |data$pred$mappabiliy_40mer < 1 | data$pred$mappability_24mer < 1 | data$muts$mutated == 0
  #tMuts = tMuts[!toExclude,]
  mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
  rm(tMuts, data)

  # load predictions
  preds = do.call(rbind,sapply(1:31, function(i){
     # RF predictions
     temp = get(load(paste0("data/procData/exomeMuts/exomeMuts_GLMpredictions_part",
                            i, "_", tissue, ".RData")))
     colnames(temp)[colnames(temp) == "prediction"] = "GLMpreds"
     # context predictions
     temp$contextPreds = get(load(paste0("data/procData/exomeMuts/exomeMuts_contextPredictions_part",
                                         i, "_", tissue, ".RData")))$prediction
     # get info which positions were mutated
     allPos = paste(temp$chr, temp$pos, sep="_")
     temp$mutated = allPos %in% mutPos
     return(temp[,c("mutated", "GLMpreds", "contextPreds")])
  }, simplify =F))
  save(preds, file = paste0("data/procData/exomeMuts/combinedPredictions_",
                            tissue, ".RData"))
  return(NA)
})
#####
# 
# 
# # create combinations of the two scores and compare using AUC ######
# print("testing combinations")
# png("fig/exome_predCombinations_AUC.png",
#     height = 2000,width = 1400,pointsize=30)
# par(mfrow = c(5,2), mar = c(2.5,9,1,0.5))
# combPerfs = sapply(tissues, function(tissue){
#    print(tissue)
#    load(paste0("data/procData/exomeMuts/combinedPredictions_",
#                tissue, ".RData"))
#    # preds = preds[sample(1:nrow(preds),1000000),]
# 
# 
#    # compute corrected scores
#    temp = prediction(predictions=preds$RFpreds, labels=preds$mutated)
#    temp = performance(temp,"prec")
#    temp = approxfun(x = temp@x.values[[1]], y = temp@y.values[[1]])
#    preds$RFcorr = temp(preds$RFpreds)
#    temp = prediction(predictions=preds$contextPreds, labels=preds$mutated)
#    temp = performance(temp,"prec")
#    temp = approxfun(x = temp@x.values[[1]], y = temp@y.values[[1]])
#    preds$contextCorr = temp(preds$contextPreds)
# 
#    # compute scaled RF score
#    preds$RFscaled = (preds$RFpreds/0.5)*(sum(preds$mutated)/nrow(preds))
# 
#    # add simple combinations
#    preds$mult = preds$contextPreds * preds$RFpreds
#    preds$mean = (preds$contextPreds + preds$RFpreds)/2
#    temp = preds$RFpreds/(1-preds$RFpreds) * preds$contextPreds/(1-preds$contextPreds)
#    preds$combOdds = temp/(1+temp)
# 
#    # add combinations based on rescaled RF predictions
#    preds$multScaled = preds$contextPreds * preds$RFscaled
#    preds$meanScaled = (preds$contextPreds + preds$RFscaled)/2
#    temp = preds$RFscaled/(1-preds$RFscaled) * preds$contextPreds/(1-preds$contextPreds)
#    preds$combOddsScaled = temp/(1+temp)
# 
#    # add combinations based on corrected RF predictions
#    preds$multCorr = preds$RFcorr * preds$contextCorr
#    preds$meanCorr  = (preds$RFcorr + preds$contextCorr)/2
#    temp = preds$RFcorr/(1-preds$RFcorr) * preds$contextCorr/(1-preds$contextCorr)
#    preds$combOddsCorr= temp/(1+temp)
# 
# 
#    rm(temp); gc()
#    # compute ROCs
#    aucs = sapply(c("RFpreds", "RFscaled", "RFcorr",
#                    "contextPreds","contextCorr",
#                    "mult", "mean", "combOdds",
#                    "multScaled", "meanScaled", "combOddsScaled",
#                    "multCorr", "meanCorr", "combOddsCorr"),
#                  function(x){
#                     pred = prediction(predictions = preds[,x],
#                                       labels=preds$mutated)
#                     performance(pred, "auc")@y.values[[1]]
#                  })
#    cols =  setNames(hcl.colors(n=length(aucs), palette = "viridis"),names(sort(aucs)))
#    barplot(aucs, horiz = T, las = 1, col = cols[names(aucs)], main = tissue, xlab = "AUC")
#    return(aucs)
# }, simplify = F)
# dev.off()
# #####
# # --> multCorr is the best across tissues, combOdds is the best without overtraining danger.
# 
#  
# # prepare tables with chr, pos, RFpred, contextPred, combOdds, combOddsCorr ######
# print("prepare combinations")
# dumpVar = sapply(tissues, function(tissue){
#    print(tissue)
#    # load information which positions were mutated
#    load(paste0("data/procData/traindata/traindata_", tissue, ".RData"))
#    tMuts = data$muts
#    toExclude = data$pred$ConsensusExcludable == 1 |
#       data$pred$repeatMasker ==1 |
#       data$pred$tandemRepeatFinder == 1 |
#       data$muts$mutated == 0
#    tMuts = tMuts[!toExclude,]
#    mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
#    rm(tMuts, data, toExclude)
#    # load predictions
#    preds = do.call(rbind,sapply(1:32, function(i){
#       # RF predictions
#       temp = get(load(paste0("data/procData/exomeMuts/exomeMuts_RFpredictions_part",
#                              i, "_", tissue, ".RData")))
#       colnames(temp)[colnames(temp) == "prediction"] = "RFpreds"
#       # context predictions
#       temp$contextPreds = get(load(paste0("data/procData/exomeMuts/exomeMuts_contextPredictions_part",
#                                           i, "_", tissue, ".RData")))$prediction
#       allPos = paste(temp$chr, temp$pos, sep="_")
#       temp$mutated = allPos %in% mutPos
#       return(temp[,c("chr", "pos", "mutated", "RFpreds", "contextPreds")])
#    }, simplify =F))
#    # create combinations combOdds and combOddsCorr
#    temp = prediction(predictions=preds$RFpreds, labels=preds$mutated)
#    temp = performance(temp,"prec")
#    temp = approxfun(x = temp@x.values[[1]], y = temp@y.values[[1]])
#    preds$RFcorr = temp(preds$RFpreds)
#    temp = prediction(predictions=preds$contextPreds, labels=preds$mutated)
#    temp = performance(temp,"prec")
#    temp = approxfun(x = temp@x.values[[1]], y = temp@y.values[[1]])
#    preds$contextCorr = temp(preds$contextPreds)
#    temp = preds$RFpreds/(1-preds$RFpreds) * preds$contextPreds/(1-preds$contextPreds)
#    preds$combOdds = temp/(1+temp)
#    temp = preds$RFcorr/(1-preds$RFcorr) * preds$contextCorr/(1-preds$contextCorr)
#    preds$combOddsCorr = temp/(1+temp)
# 
#    save(preds, file = paste0("data/procData/exomeMuts/exomePredictions_",
#                              tissue, ".RData"))
# })
# ######
# 
# # compare PRs and ROCs #####
# print("creating PRs and ROCs")
# cols = setNames(rainbow(6),
#                 c("RFpreds", "contextPreds", 
#                   "RFcorr", "contextCorr", "combOdds", "combOddsCorr"))
# png("fig/exome_predCombinations_PR.png", 
#     height = 2000,width = 2000,pointsize=30)
# par(mfrow = c(5,4), mar = c(3,3,1,0.5))
# dumpVar = sapply(tissues, function(tissue){
#    print(tissue)
#    load(paste0("data/procData/exomeMuts/exomePredictions_",
#                tissue, ".RData"))
#    # preds = preds[sample(1:nrow(preds), 1000000),]
#    
#    pred = prediction(predictions = preds[,-(1:3)],
#                      labels=list(preds$mutated,preds$mutated,preds$mutated,
#                                  preds$mutated,preds$mutated,preds$mutated))
#    # plot ROC
#    perf = performance(pred, "tpr", "fpr")
#    plot(perf, col = as.list(cols[colnames(preds)[-(1:3)]]),
#         lwd = 2, main = tissue, mgp = c(2,0.5,0))
#    abline(0,1, col = "grey", lty = 2)
#    legend("bottomright", col = cols[colnames(preds)[-(1:3)]], 
#           legend = names(cols[colnames(preds)[-(1:3)]]), lwd = 2,
#           cex = 0.7)
#    # plot PR
#    perf = performance(pred, "prec", "rec")
#    plot(perf, col = as.list(cols[colnames(preds)[-(1:3)]]),
#         lwd = 2, main = tissue, mgp = c(2,0.5,0))
#    legend("topright", col = cols[colnames(preds)[-(1:3)]], 
#           legend = names(cols[colnames(preds)[-(1:3)]]), lwd = 2,
#           cex = 0.7)
#    return(NA)
# }, simplify = F)
# dev.off()
# #####


# GTEx mutations ######
print("GTEx")
cols = setNames(brewer.pal(n = 4, name = "Set1"),
                c("RFpreds", "contextPreds",
                  "combOddsCorr", "combOdds"))
png("fig/exome_GTExPerformance_ROCPR.png",
    height = 2000,width = 2000,pointsize=30)
par(mfrow = c(5,4), mar = c(3,3,1,0.5))
combPerfs = sapply(tissues, function(tissue){
   print(tissue)
   # load predictions
   load(paste0("data/procData/exomeMuts/exomePredictions_",
               tissue, ".RData")) #preds

   # load information which positions were mutated
   load(paste0("data/procData/validationData/validationData_GTEx_", tissue, ".RData"))
   tMuts = data$muts
   toExclude = data$pred$ConsensusExcludable == 1 |
      data$pred$repeatMasker ==1 |
      data$pred$tandemRepeatFinder == 1 |
      data$muts$mutated == 0
   tMuts = tMuts[!toExclude,]
   mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
   rm(tMuts, data, toExclude)

   # get info which positions were mutated
   allPos = paste(preds$chr, preds$pos, sep="_")
   preds$mutated = allPos %in% mutPos
   rm(allPos, mutPos); gc()

   # get performance
   preds = preds[c(which(preds$mutated), sample(which(!preds$mutated),  100000-sum(preds$mutated))),]
   pred = prediction(predictions = preds[,c("RFpreds", "contextPreds",
                                            "combOddsCorr", "combOdds")],
                     labels=list(preds$mutated,preds$mutated,
                                 preds$mutated,preds$mutated))

   # plot ROC
   perf = performance(pred, "tpr", "fpr")
   plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0),
        col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]))
   legend("bottomright", lwd = 2,
          legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
          col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
   abline(0,1,lty = 2, col = "grey", lwd = 2)

   # plot PR
   perf = performance(pred, "prec", "rec")
   plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0),
        col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]),)
   legend("topright", lwd = 2,
          col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
          legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"))

   aucs  = performance(pred,"auc")
   aucs = unlist(aucs@y.values)
   names(aucs) = c("RFpreds", "contextPreds",
                   "combOddsCorr", "combOdds")
   return(aucs)
})
dev.off()
save(combPerfs, file = "data/rdata/exome_GTExPerformance_AUC.RData")
# plot AUCs
png("fig/exome_GTExPerformance_AUC.png",
    height = 1000,width = 2000,pointsize=30)
par(mar = c(3,4,1.5,8))
barplot(combPerfs, beside = T, ylab = "AUC",las = 1, mgp = c(2.5,1,0),
        col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
        main = "GTEx mutations")
legend(51,0.8, xpd=NA,
       legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
       fill =  cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
dev.off()
#####


# single tissue datasets #####
print("single tissue datasets")
png("fig/exome_HealthyTissuePerformance_ROCPR.png", 
    height = 1700,width = 2000,pointsize=30)
par(mfrow = c(4,4), mar = c(3,3,1,0.5))
combPerfs = sapply(
   c("brain", "colon","esophagus", "liver","liver_Blokzijl",
     "lung", "prostate", "skin"), 
   function(tissue){
      print(tissue)
      # load predictions
      if(tissue == "liver_Blokzijl"){
         load(paste0("data/procData/exomeMuts/exomePredictions_",
                     "liver", ".RData")) 
      } else{
         load(paste0("data/procData/exomeMuts/exomePredictions_",
                     tissue, ".RData")) #preds
      }
      
      # preds = preds[sample(1:nrow(preds), 1000000),]
      
      # load information which positions were mutated
      load(paste0("data/procData/validationData/validationData_",
                  tissue, ".RData"))
      tMuts = data$muts
      toExclude = data$pred$ConsensusExcludable == 1 | 
         data$pred$repeatMasker ==1 | 
         data$pred$tandemRepeatFinder == 1 | 
         data$muts$mutated == 0
      tMuts = tMuts[!toExclude,]
      mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
      rm(tMuts, data, toExclude)
      
      # get info which positions were mutated
      allPos = paste(preds$chr, preds$pos, sep="_")
      preds$mutated = allPos %in% mutPos
      
      pred = prediction(predictions = preds[,c("GLMpreds", "contextPreds",
                                               "combOddsCorr", "combOdds")],
                        labels=list(preds$mutated,preds$mutated,
                                    preds$mutated,preds$mutated))
      
      # plot ROC
      perf = performance(pred, "tpr", "fpr")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0),
           col = as.list(cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")]))
      legend("bottomright", lwd = 2,
             legend = c("GLMpreds", "contextPreds","combOddsCorr", "combOdds"),
             col = cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")])
      abline(0,1,lty = 2, col = "grey", lwd = 2)
      
      # plot PR
      perf = performance(pred, "prec", "rec")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0), 
           col = as.list(cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")]),)
      legend("topright", lwd = 2, 
             col = cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")],
             legend = c("GLMpreds", "contextPreds","combOddsCorr", "combOdds"))
      
      aucs  = performance(pred,"auc")
      aucs = unlist(aucs@y.values)
      names(aucs) = c("GLMpreds", "contextPreds",
                      "combOddsCorr", "combOdds")
      return(aucs)
   })
dev.off()
save(combPerfs, file = "data/rdata/exome_HealthyTissuePerformance_AUC.RData")

# plot AUCs

png("fig/exome_HealthyTissuePerformance_AUC.png", 
    height = 1000,width = 2000,pointsize=30)
par(mar = c(3,4,1.5,8))
barplot(combPerfs, beside = T, ylab = "AUC",las = 1, mgp = c(2.5,1,0),
        col = cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")],
        main = "healthy tissue mutations, mixed sources")
legend(41,0.6, xpd=NA,
       legend = c("GLMpreds", "contextPreds","combOddsCorr", "combOdds"),
       fill =  cols[c("GLMpreds", "contextPreds","combOddsCorr", "combOdds")])
dev.off()



#####

# Moore et al. mutations ######
print("Moore")
png("fig/exome_MoorePerformance_ROCPR.png", 
    height = 1700,width = 2000,pointsize=30)
par(mfrow = c(4,4), mar = c(3,3,1,0.5))
combPerfs = sapply(
   c("colon","esophagus", "kidney","liver", "lung", 
     "prostate", "skin"), 
   function(tissue){
      print(tissue)
      # load predictions
      load(paste0("data/procData/exomeMuts/exomePredictions_",
                  tissue, ".RData")) #preds
      # preds = preds[sample(1:nrow(preds), 1000000),]
      
      # load information which positions were mutated
      load(paste0("data/procData/validationData/validationData_Moore_", tissue, ".RData"))
      tMuts = data$muts
      toExclude = data$pred$ConsensusExcludable == 1 | 
         data$pred$repeatMasker ==1 | 
         data$pred$tandemRepeatFinder == 1 | 
         data$muts$mutated == 0
      tMuts = tMuts[!toExclude,]
      mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
      rm(tMuts, data, toExclude)
      
      # get info which positions were mutated
      allPos = paste(preds$chr, preds$pos, sep="_")
      preds$mutated = allPos %in% mutPos
      
      pred = prediction(predictions = preds[,c("RFpreds", "contextPreds",
                                               "combOddsCorr", "combOdds")],
                        labels=list(preds$mutated,preds$mutated,
                                    preds$mutated,preds$mutated))
      
      # plot ROC
      perf = performance(pred, "tpr", "fpr")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0),
           col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]))
      legend("bottomright", lwd = 2, 
             legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
             col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
      abline(0,1,lty = 2, col = "grey", lwd = 2)
      
      # plot PR
      perf = performance(pred, "prec", "rec")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0), 
           col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]),)
      legend("topright", lwd = 2, 
             col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
             legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"))
      
      aucs  = performance(pred,"auc")
      aucs = unlist(aucs@y.values)
      names(aucs) = c("RFpreds", "contextPreds",
                      "combOddsCorr", "combOdds")
      return(aucs)
   })
dev.off()
save(combPerfs, file = "data/rdata/exome_MoorePerformance_AUC.RData")

# plot AUCs
png("fig/exome_MoorePerformance_AUC.png", 
    height = 1000,width = 2000,pointsize=30)
par(mar = c(3,4,1.5,8))
barplot(combPerfs, beside = T, ylab = "AUC",las = 1, mgp = c(2.5,1,0),
        col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
        main = "healthy tissue mutations, Moore et al.")
legend(36,0.6, xpd=NA,
       legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
       fill =  cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
dev.off()


#####

# Li et al. mutations #####
print("Li")
png("fig/exome_LiPerformance_ROCPR.png", 
    height = 1200,width = 2000,pointsize=30)
par(mfrow = c(2,4), mar = c(3,3,1,0.5))
combPerfs = sapply(
   c("colon","esophagus", "liver", "lung"), 
   function(tissue){
      print(tissue)
      # load predictions
      load(paste0("data/procData/exomeMuts/exomePredictions_",
                  tissue, ".RData")) #preds
      # preds = preds[sample(1:nrow(preds), 1000000),]
      
      # load information which positions were mutated
      load(paste0("data/procData/validationData/validationData_Li_", tissue, ".RData"))
      tMuts = data$muts
      toExclude = data$pred$ConsensusExcludable == 1 | 
         data$pred$repeatMasker ==1 | 
         data$pred$tandemRepeatFinder == 1 | 
         data$muts$mutated == 0
      tMuts = tMuts[!toExclude,]
      mutPos = paste(tMuts$chr, tMuts$pos, sep="_")
      rm(tMuts, data, toExclude)
      
      # get info which positions were mutated
      allPos = paste(preds$chr, preds$pos, sep="_")
      preds$mutated = allPos %in% mutPos
      
      pred = prediction(predictions = preds[,c("RFpreds", "contextPreds",
                                               "combOddsCorr", "combOdds")],
                        labels=list(preds$mutated,preds$mutated,
                                    preds$mutated,preds$mutated))
      
      # plot ROC
      perf = performance(pred, "tpr", "fpr")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0),
           col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]))
      legend("bottomright", lwd = 2,
             legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
             col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
      abline(0,1,lty = 2, col = "grey", lwd = 2)
      
      # plot PR
      perf = performance(pred, "prec", "rec")
      plot(perf, lwd = 2, main = tissue, mgp = c(2,0.5,0), 
           col = as.list(cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")]),)
      legend("topright", lwd = 2, 
             col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
             legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"))
      
      aucs  = performance(pred,"auc")
      aucs = unlist(aucs@y.values)
      names(aucs) = c("RFpreds", "contextPreds",
                      "combOddsCorr", "combOdds")
      return(aucs)
   })
dev.off()
save(combPerfs, file = "data/rdata/exome_LiPerformance_AUC.RData")

# plot AUCs
png("fig/exome_LiPerformance_AUC.png", 
    height = 1000,width = 2000,pointsize=30)
par(mar = c(3,4,1.5,8))
barplot(combPerfs, beside = T, ylab = "AUC",las = 1, mgp = c(2.5,1,0),
        col = cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")],
        main = "healthy tissue mutations, Li et al.")
legend(21,0.7, xpd=NA,
       legend = c("RFpreds", "contextPreds","combOddsCorr", "combOdds"),
       fill =  cols[c("RFpreds", "contextPreds","combOddsCorr", "combOdds")])
dev.off()



######


dumpVar = sapply(c("GTEx", "VariousTissues", "Moore", "Li"), function(datsource){
   print(datsource)
   if(datsource == "GTEx"){
      load("data/rdata/exome_GTExPerformance_AUC.RData")
      load("data/rdata/GTEx_nPerSample.RData")
   } else if(datsource == "VariousTissues"){
      load("data/rdata/exome_HealthyTissuePerformance_AUC.RData")
      load("data/rdata/healthyTissues_nPerSample.RData")
   } else if (datsource == "Moore"){
      load("data/rdata/exome_MoorePerformance_AUC.RData")
      load("data/rdata/Moore_nPerSample.RData")
   }else if (datsource == "Li"){
      load("data/rdata/exome_LiPerformance_AUC.RData")
      load("data/rdata/Li_nPerSample.RData")
   }
   png(paste0("fig/exome_", datsource, "Performance_AUCvsNperSample.png"),
       width=800, height=500, pointsize=15)
   suppressWarnings({
      par(mfrow = c(2,3), mar = c(4,4,0.5,0.5), oma = c(0,0,3,0))
      tissueCols = setNames(rainbow(length(nPerSample)), names(nPerSample))
      AUCs = combPerfs["RFpreds",names(nPerSample)]
      # plot 1: AUC vs. total number of mutations
      nTotal = sapply(nPerSample, function(x){
         sum(sapply(x, sum))
      })
      plot(nTotal, AUCs, mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n mutations per tissue", ylab = "RF AUC", 
           bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nTotal), col = "grey", lty = 2)
      # AUC vs. n muts per sample
      nMutsPerSample = sapply(nPerSample, unlist)
      nMean = sapply(nMutsPerSample, mean)
      nMin = sapply(nMutsPerSample, min)
      nMax = sapply(nMutsPerSample, max)
      plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n mutations per sample", ylab = "RF AUC",
           bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
      arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
             y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
      points(nMean, AUCs, bg = tissueCols[names(AUCs)], pch = 21)
      # text(x = nMean, y = AUCs, labels=names(AUCs), 
      #      col = tissueCols, pos=1)
      # AUC vs. n samples per tissue
      nSampTotal = sapply(nPerSample, function(x){
         sum(sapply(x, length))
      })
      plot(nSampTotal, AUCs, mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n samples per tissue", ylab = "RF AUC",
           bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nSampTotal), col = "grey", lty = 2)
      # AUC vs. n Samples per Donor
      nSamp = sapply(nPerSample, function(x){
         sapply(x, length)
      })
      nMean = sapply(nSamp, mean)
      nMin = sapply(nSamp, min)
      nMax = sapply(nSamp, max)
      plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n samples per donor", ylab = "RF AUC", 
           bg = tissueCols[names(AUCs)], pch = 21)
      arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
             y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
      points(nMean, AUCs, bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
      # AUC vs. n muts per donor
      nMutsPerDonor = sapply(nPerSample, function(x){sapply(x, sum)})
      nMean = sapply(nMutsPerDonor, mean)
      nMin = sapply(nMutsPerDonor, min)
      nMax = sapply(nMutsPerDonor, max)
      plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n mutations per donor", ylab = "RF AUC",
           bg = tissueCols[names(AUCs)], pch = 21)
      arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
             y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
      points(nMean, AUCs, bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
      # AUC vs. n donor per tissue
      nDonor = sapply(nPerSample, length)
      plot(nDonor, AUCs, mgp = c(2.5,1,0),
           ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
           xlab = "n donors per tissue", ylab = "RF AUC",  
           bg = tissueCols[names(AUCs)], pch = 21)
      abline(lm(AUCs ~ nDonor), col = "grey", lty = 2)
      
      legend("topright", pt.bg = tissueCols[names(AUCs)], pch = 21, 
             legend=names(AUCs), ncol = 2)
      title(main = datsource, outer = T)
   })
   dev.off()
})


png(paste0("fig/exome_MooreAndLi_Performance_AUCvsNperSample.png"),
    width=800, height=500, pointsize=15)
MoorePerf = get(load("data/rdata/exome_MoorePerformance_AUC.RData"))
MooreNsamp = get(load("data/rdata/Moore_nPerSample.RData"))
LiPerf = get(load("data/rdata/exome_LiPerformance_AUC.RData"))
LiNsamp = get(load("data/rdata/Li_nPerSample.RData"))

par(mfrow = c(2,3), mar = c(4,4,0.5,0.5), oma = c(0,0,3,0))
tissueCols = setNames(rainbow(length(tissues)), tissues)
AUCs = c(MoorePerf["RFpreds",names(MooreNsamp)],
         LiPerf["RFpreds",names(LiNsamp)])
source = c(rep("Moore", ncol(MoorePerf)), rep("Li", ncol(LiPerf)))
# plot 1: AUC vs. total number of mutations
nTotal = c(sapply(MooreNsamp, function(x){sum(sapply(x, sum))}), 
           sapply(LiNsamp, function(x){sum(sapply(x, sum))}))
plot(nTotal, AUCs, mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n mutations per tissue", ylab = "RF AUC", cex = 1.5,
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
abline(lm(AUCs ~ nTotal), col = "grey", lty = 2)
# AUC vs. n muts per sample
nMutsPerSample = c(sapply(MooreNsamp, unlist),sapply(LiNsamp, unlist))
nMean = sapply(nMutsPerSample, mean)
nMin = sapply(nMutsPerSample, min)
nMax = sapply(nMutsPerSample, max)
plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n mutations per sample", ylab = "RF AUC", cex = 1.5,
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
       y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
points(nMean, AUCs, bg = tissueCols[names(AUCs)],
       pch = c(21,24)[(source=="Li") + 1], cex = 1.5 )
# AUC vs. n samples per tissue
nSampTotal = c(sapply(MooreNsamp, function(x){sum(sapply(x, length))}),
              sapply(LiNsamp, function(x){sum(sapply(x, length))}))
plot(nSampTotal, AUCs, mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n samples per tissue", ylab = "RF AUC", cex = 1.5,
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
abline(lm(AUCs ~ nSampTotal), col = "grey", lty = 2)
# AUC vs. n Samples per Donor
nSamp = c(sapply(MooreNsamp, function(x){sapply(x, length)}),
         sapply(LiNsamp, function(x){sapply(x, length)}))
nMean = sapply(nSamp, mean)
nMin = sapply(nSamp, min)
nMax = sapply(nSamp, max)
plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n samples per donor", ylab = "RF AUC",  cex = 1.5,
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
       y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
points(nMean, AUCs, bg = tissueCols[names(AUCs)],
       pch = c(21,24)[(source=="Li") + 1], cex = 1.5)
abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
# AUC vs. n muts per donor
nMutsPerDonor = c(sapply(MooreNsamp, function(x){sapply(x, sum)}),
                  sapply(LiNsamp, function(x){sapply(x, sum)}))
nMean = sapply(nMutsPerDonor, mean)
nMin = sapply(nMutsPerDonor, min)
nMax = sapply(nMutsPerDonor, max)
plot(nMean, AUCs, xlim = c(min(nMin), max(nMax)),mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n mutations per donor", ylab = "RF AUC", cex = 1.5,
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
arrows(x0 = nMin[nMax-nMin != 0], x1 = nMax[nMax-nMin != 0],
       y0=AUCs[nMax-nMin != 0], angle = 90, length = 0.03, code = 3)
points(nMean, AUCs, bg = tissueCols[names(AUCs)], 
       pch = c(21,24)[(source=="Li") + 1], cex = 1.5)
abline(lm(AUCs ~ nMean), col = "grey", lty = 2)
# AUC vs. n donor per tissue
nDonor = c(sapply(MooreNsamp, length), sapply(LiNsamp, length))
plot(nDonor, AUCs, mgp = c(2.5,1,0),
     ylim = c(min(AUCs)-((max(AUCs)-min(AUCs))*0.1), max(AUCs)), las = 1,
     xlab = "n donors per tissue", ylab = "RF AUC",  cex = 1.5, 
     bg = tissueCols[names(AUCs)], pch = c(21,24)[(source=="Li") + 1])
abline(lm(AUCs ~ nDonor), col = "grey", lty = 2)

legend("topright", pt.bg = tissueCols, pch = 21, 
       legend=names(tissueCols), ncol = 2, pt.cex =1.5)
legend("bottomright", pt.bg = "grey", pch = c(21,24), legend=c("Moore", "Li"), 
       pt.cex = 1.5)
title(main = "Moore and Li data sizes", outer = T)
dev.off()
