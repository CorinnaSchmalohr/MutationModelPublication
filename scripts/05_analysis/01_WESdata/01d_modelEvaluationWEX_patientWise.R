
# preparation #####
library(ROCR)
source("scripts/05_analysis/00_NamesAndColors.R")
plotEnding = "_20230916"
source("lib/general_function.R")

bases = c("A", "C", "G", "T")
mutTypes = cbind(paste0(paste0(rep(rep(bases, each = 4), 6),"[",
                                 c(rep("C", 48), rep("T", 48)), ">",
                                 c(rep(c("A", "G", "T"), each =16),
                                   rep(c("A", "C", "G"), each =16)), "]",
                                 rep(bases, 24))),
                   paste0(paste0(rep(rep(rev(bases), each = 4),  6),"[",
                                 c(rep("G", 48), rep("A", 48)), ">",
                                 c(rep(c("T", "C", "A"), each =16),
                                   rep(c("T", "G", "C"), each =16)), "]",
                                 rep(rev(bases), 24))))
rownames(mutTypes) = paste(mutTypes[,1])
mutTypeTranslator= setNames(nm=c(mutTypes[,1], mutTypes[,2]),
                              object=c(rownames(mutTypes), rownames(mutTypes)))
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
cancer2Tissue = do.call(rbind,sapply(names(tissue2Cancer), function(tissue){
  cbind(tissue, tissue2Cancer[[tissue]])
}))
cancer2Tissue = setNames(object = cancer2Tissue[,1], nm = cancer2Tissue[,2])
chrs = paste0("chr", c(1:22))
signatures = read.table("data/rawdata/COSMIC_catalogue-signatures_SBS96_v3.4/COSMIC_v3.4_SBS_GRCh37.txt", 
                        header = T, row.names = 1)
signatures = signatures[mutTypes[,1],]
######


# get table with relevant mutation  data #####
load("data/MutTables/exomeTrainData/muts.RData")
mutData = muts[,c("Chromosome", "Start_Position","Reference_Allele", "Tumor_Seq_Allele2", "Tumor_Sample_Barcode", "CONTEXT","cancerType")]
mutData$patientBarcodes = substr(mutData$Tumor_Sample_Barcode, 1,12)
mutData$mutType = paste0(substr(mutData$CONTEXT,5,5),
                      "[",mutData$Reference_Allele, ">",mutData$Tumor_Seq_Allele2, "]",
                      substr(mutData$CONTEXT,7,7))
mutData$tissue = cancer2Tissue[mutData$cancerType]
mutData$position = paste0(mutData$Chromosome, "_",mutData$Start_Position)
#####


# get patient information #####
meta = read.table("~/Documents/clinical_PANCAN_patient_with_followup.tsv",
                  header = T, sep = "\t", quote = "")
table(mutData$patientBarcodes %in% meta$bcr_patient_barcode) # not all patients available, no idea why
toRemove = sapply(meta, function(x){
  all(x =="" | x == "[Not Applicable]" | x == "[Not Available]", na.rm = T) })
meta = meta[,!toRemove]

subMeta = meta[meta$bcr_patient_barcode %in% mutData$patientBarcodes,c("bcr_patient_barcode", "days_to_birth", "radiation_therapy")]
rownames(subMeta) = subMeta[,1]
#####


# count mutations per sample  #####
png(paste0("fig/nMutsPerSampleDistribution", plotEnding, ".png"), 
    width = 1200, height = 700, pointsize = 15)
par(mfrow = c(2,5))
topPercent = sapply(tissues, function(tissue){
  temp = table(mutData$patientBarcodes[mutData$tissue == tissue])
  temp = sort(temp, decreasing = T)
  barplot(temp,main = tissue, las = 1, ylab = "n Mutations", xaxt = "n", xlab = "samples")
  topPerc = temp[1:5]/sum(temp)*100
  legend("topright", format(topPerc, digits = 2), title = "percent of top 5")
  topPerc
}, simplify = F)
dev.off()

top5Table = sapply(tissues, function(tissue){
  temp = table(mutData$patientBarcodes[mutData$tissue == tissue])
  temp = sort(temp, decreasing = T)
  temp = temp/sum(temp)
  c(temp[1:5], "Rest" = sum(temp[-(1:5)]))
})
top5Table = top5Table[,rev(tissues)]
pdfAndPng(file = paste0("fig/nMutsPerSampleCompositions", plotEnding), 
          width = 12, height = 5,pngArgs = list(pointsize=20), 
          pdfArgs = list(pointsize=15), expr = expression({
            par(mar = c(3.5,6,1,8))
            temp = barplot(top5Table, las = 1, horiz = T, col = c(hcl.colors(5, palette = "Dark 3"), "grey"), xlab = "Proportion of Mutations",
                           mgp = c(2.5,1,0),
                           names.arg = t2T[colnames(top5Table)], density = c(rep(NA,5), 20))
            legend(x = 1, y = ncol(top5Table)+2, xjust = -0.1, yjust = 1,
                   legend = c(paste0("Sample ", 1:5), "Rest"), 
                   fill = c(2:6, "grey"),xpd = T, density = c(rep(NA,5), 40))
            toMark = which(top5Table[1:5,] >= 0.05, arr.ind = T)
            dumpVar = apply(toMark, 1,function(indic){
              tis = indic[2]
              pat = indic[1]
              val = round(top5Table[pat,tis]*100, digits = 1)
              xCoords = c(0,cumsum(top5Table[,tis]))
              text(x = mean(xCoords[c(pat, pat+1)]), y = temp[tis], cex =0.8,
                   labels = paste0(val, "%"))
            })
}))




outliers = do.call(rbind,sapply(tissues, function(tissue){
  temp = topPercent[[tissue]]
  res = temp[temp>=5]
  if(length(res)>0){
    return(cbind(tissue, names(res), res))
  } else{
    NULL
  }
}))
# meta[meta$bcr_patient_barcode %in% outliers[,2],] --> nothing interesting
outlierInfo  = cbind(outliers, subMeta[outliers[,2],])
outlierInfo
#####


# count mutation types per sample #####
dumpVar = sapply(unique(outliers[,"tissue"]), function(tissue){
  subDat = mutData[mutData$tissue == tissue,]
  # first for entire dataset
  subSig = factor(mutTypeTranslator[subDat$mutType], levels = mutTypes[,1])
  png(paste0("fig/mutationOutliers_", tissue, ".png"), width = 1000, height = 1000, pointsize = 15)
  par(mfrow = c(1,sum(outliers[,1] == tissue)*3+1), mar = c(2,3,2,0), oma = c(2,4,0,0))
  barplotTemp = barplot(table(subSig), las = 1, horiz = T, cex.names = 0.6,
                        main = "All samples", names.arg = NA, col = "black")
  axis(side = 2, at = barplotTemp, labels = names(table(subSig)),
       family = "mono", las = 1, cex = 0.8)
  # then for each of the two most prominent patients
  patients = outliers[outliers[,1] == tissue,2]
  patientMutCounts = sapply(patients, function(patient){
    patientMuts = factor(mutTypeTranslator[subDat$mutType[subDat$patientBarcodes == patient]], 
                         levels = mutTypes[,1])
    patientCol = which(patients == patient)+1
    barplot(table(patientMuts), las = 1, horiz = T, yaxt = "n", main = patient, 
            col = patientCol)
    mostCorrSign = cor(table(patientMuts), signatures)
    maxCorr = which.max(mostCorrSign)
    barplot(signatures[,maxCorr], las = 1, horiz = T, yaxt = "n", 
            main = paste0("\n",colnames(signatures)[maxCorr], "\nr = ", 
                          format(mostCorrSign[maxCorr], digits = 2)))
    nonPatientMuts = factor(mutTypeTranslator[subDat$mutType[subDat$patientBarcodes != patient]], 
                            levels = mutTypes[,1])
    barplot(table(nonPatientMuts), las = 1, horiz = T, yaxt = "n",
            main = paste0("\nWithout\n", patient), col = addAlpha(patientCol, alpha = 0.5))
    return(table(patientMuts))
  })
  # patientMutCounts = cbind("Mutation Types" = rownames(patientMutCounts), patientMutCounts)
  # write.table(patientMutCounts,
  #             file = paste0("data/MutTables/exomeTrainData/",tissue, "_outliers.txt"), 
  #             quote = F, row.names = F, sep = "\t")
  mtext("Count", side = 1, outer = T, line = 1)
  mtext("Mutation type", side = 2, outer = T, line = 2.5)
  dev.off()
})
#####

# Compute prediction performance stratified by (outlier) patients #####
# here we have predictions for each tissue (per chromosome)
load("data/Modeling/exomeTrainData/RF/predPerTissueRF.RData")
# we do not have the information of the patient in our processed data, therefore
# we have to infer them based on chromosomal position
load("data/processedData/pentamers.RData")

perfPatients = sapply(tissues, function(tissue){
  print(tissue)
  # load the preprocessed data where we have the positions in same order as predictions
  load(paste0("data/MutTables/exomeTrainData/", tissue, "_Muts_mapped.RData")) #data
  # translation from position to patient of origin (only possible for TPs of course)
  subMutData = mutData[mutData$tissue == tissue,c("patientBarcodes", "position")]
  multiple = duplicated(subMutData$position) | duplicated(subMutData$position, fromLast=T)
  subMutData = subMutData[!multiple,]
  pos2Patient = setNames(object = subMutData$patientBarcodes,subMutData$position)
  # combine mutation position with predictions
  subPreds = cbind(data$muts,do.call(rbind,predPerTissueRF[[tissue]]))
  # now use position to infer patient
  subPreds$mutation = paste0(subPreds$chr, "_", subPreds$pos)
  subPreds$patient = pos2Patient[subPreds$mutation]
  
  # since TNs don't come from a patient, split them across patients according to sequence context
  for(context in rownames(fivemers)){
    subPreds$patient[(subPreds$context %in% fivemers[context,]) & subPreds$label == 0] = 
      sample(subPreds$patient[(subPreds$context %in% fivemers[context,]) & subPreds$label == 1])
  }
  patientCounts = table(subPreds$patient)
  mostFrequentPatients = names(sort(patientCounts, decreasing = T))[1:5]
  perfs = sapply(mostFrequentPatients, function(patient){
    x = subPreds[subPreds$patient != patient,]
    perf = prediction(x$pred, x$label)
    roc = performance(perf, "tpr", "fpr")
    pr = performance(perf,"prec", "rec")
    auc = performance(perf,"auc")@y.values[[1]]
    return(list(roc = roc, pr = pr, auc = auc))
  }, simplify = F)
  return(perfs)
}, simplify = F)
#####


# visualize prediction performance  stratified by (outlier) patients ######
# load performances on whole data
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
# AUROC
AUCpatients = sapply(tissues, function(tissue){
  x = perfPatients[[tissue]]
  c(wholeData = ROC_PR_RF_concat[[tissue]]$auc@y.values[[1]], sapply(x, function(y){y$auc}))
})
rownames(AUCpatients) = c("Whole data", paste0("Without Sample ", 1:5))
pdfAndPng(file = paste0("fig/nMutOutlierPerformanceAUROC", plotEnding), 
          width = 9, height = 5, pngArgs = list(pointsize = 15), pdfArgs  = list(pointsize = 10),
          expr = expression({
            par(mar =c(3,4,0.5,0.5))
            barplot(AUCpatients-0.5, beside = T, legend.text = T, las = 1, 
                    names.arg = t2T[tissues], ylab = "AUC", yaxt = "n", 
                    args.legend = list(x="topleft", inset = 0.02, cex = 0.8), 
                    col = c("black", hcl.colors(5, palette = "Dark 3")))
            temp = axTicks(2)
            axis(2, las = 1, at = axTicks(2), labels = c(0,axTicks(2)[-1]+0.5))
            plotrix::axis.break(axis=2, breakpos = 0.01, style = "gap")
            abline(h=0)
          }))

# dev.off()
# there is a problem for brain
#####



#  create dataset without extreme hypermutators for brain######
load(paste0("data/MutTables/exomeTrainData/brain_Muts_mapped.RData")) #data

# translation from position to patient of origin (only possible for TPs of course)
subMutData = mutData[mutData$tissue == "brain",c("patientBarcodes", "position")]
multiple = duplicated(subMutData$position) | duplicated(subMutData$position, fromLast=T)
subMutData = subMutData[!multiple,]
pos2Patient = setNames(object = subMutData$patientBarcodes,subMutData$position)
# now use position to infer patient
data$muts$mutation = paste0(data$muts$chr, "_", data$muts$pos)
data$muts$patient = pos2Patient[data$muts$mutation]
# since TNs don't come from a patient, split them across patients according to sequence context
for(context in rownames(fivemers)){
  data$muts$patient[(data$muts$context %in% fivemers[context,]) & data$muts$mutated == 0] = 
    sample(data$muts$patient[(data$muts$context %in% fivemers[context,]) & data$muts$mutated == 1])
}
datWoutHyper = cbind(data$pred[data$muts$patient != "TCGA-DU-6392",], 
            mutated = as.factor(data$muts$mutated[data$muts$patient != "TCGA-DU-6392"]))
datWoutHyperchroms = data$muts$chr[data$muts$patient != "TCGA-DU-6392"]
datWoutHyper = as.data.frame(datWoutHyper)
#####

# Retrain model with CWCV without hypermutator #####
chroms = unique(datWoutHyperchroms)
nThreads=10
library(ranger)
predictions = lapply(chroms, function(cr){
  cat(cr, ' ')
  trainData = datWoutHyper[datWoutHyperchroms != cr,]
  rf = ranger(mutated ~ ., data = trainData,
              write.forest = T, seed = 1234, num.threads =  nThreads,
              respect.unordered.factors = 'partition',
              probability = T, verbose=F, importance = "permutation")
  save(rf, file = paste0("data/Modeling/exomeTrainData/RF/", "brain", "_", 
                         cr, "_forPrediction_woutHypermutator.RData"))
  testData = datWoutHyper[datWoutHyperchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
names(predictions) = chroms
save(predictions, file = paste0("data/Modeling/exomeTrainData/RF/", "brain",
                                "_predictions_woutHypermutator.RData"))
cat('\n')


# compute CWCV performance #####
load(paste0("data/Modeling/exomeTrainData/RF/", "brain",
            "_predictions_woutHypermutator.RData"))
predConcat = do.call(rbind,predictions)
ROC = performance(prediction(predConcat$pred, 
                             predConcat$label), 
                  "tpr", "fpr")
AUC = performance(prediction(predConcat$pred, 
                             predConcat$label), 
                  "auc")
PR = performance(prediction(predConcat$pred, 
                            predConcat$label), 
                 "prec", "rec")
#####


# apply it to the full data and compute performance #####
load(paste0("data/MutTables/exomeTrainData/brain_Muts_mapped_processed.RData")) #data
predictionsFull = lapply(chroms, function(cr){
  cat(cr, ' ')
  load(paste0("data/Modeling/exomeTrainData/RF/", "brain", "_", 
                         cr, "_forPrediction_woutHypermutator.RData"))
  testData = dat[datchroms == cr,]
  p = predict(rf, data = testData, num.threads=nThreads, verbose=F)
  temp = data.frame(pred = p$predictions[,2],  label = testData$mutated)
  return(temp)
})
cat("\n")
names(predictionsFull) = chroms
predFullConcat = do.call(rbind,predictionsFull)
ROCfull = performance(prediction(predFullConcat$pred, 
                             predFullConcat$label), 
                  "tpr", "fpr")
AUCfull = performance(prediction(predFullConcat$pred, 
                             predFullConcat$label), 
                  "auc")
PRfull = performance(prediction(predFullConcat$pred, 
                            predFullConcat$label), 
                 "prec", "rec")
#####

# get the original performance ####
load("data/Modeling/exomeTrainData/RF/ROC_PR_RF_concat.RData")
ROCorig = ROC_PR_RF_concat[["brain"]]$roc
AUCorig = ROC_PR_RF_concat[["brain"]]$auc
PRorig = ROC_PR_RF_concat[["brain"]]$pr
#####


# visualize performances with and without hypermutator #####
png(paste0("fig/hypermutatorBrainRemovedPerformance", plotEnding, ".png"),
    width = 1500, height = 500, pointsize = 20)
par(mfrow = c(1,3), oma = c(0,0,0,0), mar = c(4,4,0.5,0.5))

# ROC 
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "FPR", ylab = "TPR", 
     las = 1, mgp = c(3,1,0))
abline(0,1, col = "grey", lty = 2)
plot(ROCorig, col = 1, add = T) # original
plot(perfPatients[["brain"]][["TCGA-DU-6392"]]$roc, add = T, col = 2) # wout patient in test data
plot(ROCfull, add = T, col = 3) # wout patient in train data
plot(ROC, add = T, col =  4) # wout patiernt in train and test data
legend("bottomright", col = 1:4, lty = 1, 
       legend = c("Original", "Reduced test data",
                  "Reduced train data",
                  "Reduced test and train data"))
# PR
plot(NA, xlim = c(0,1), ylim = c(0,1), xlab = "Recall", ylab = "Precision", 
     las = 1, mgp = c(3,1,0))
abline(h = 0.5, col = "grey", lty = 2)
plot(PRorig, col = 1, add = T)
plot(perfPatients[["brain"]][["TCGA-DU-6392"]]$pr, add = T, col = 2) # wout patient in test data
plot(PRfull, add = T, col = 3) # wout patient in train data
plot(PR, add = T, col =  4) # wout patiernt in train and test data
legend("bottomleft", col = 1:4, lty = 1, 
       legend = c("Original", "Reduced test data",
                  "Reduced train data",
                  "Reduced test and train data"))
# AUROC
par(mar = c(4,8,0.5,0.5))
AUCcollection = c("Original"= AUCorig@y.values[[1]],
                  "Reduced test data" = perfPatients[["brain"]][["TCGA-DU-6392"]]$auc,
                  "Reduced train data" = AUCfull@y.values[[1]],
                  "Reduced test and\ntrain data" = AUC@y.values[[1]])
barplot(rev(AUCcollection-0.5), horiz = T, las = 1, col = 4:1, xlab = "AUROC", 
        xaxt = "n", )
axis(1, at = axTicks(1), labels = axTicks(1)+0.5)
abline(v=0)
dev.off()


# AUC separately
pdfAndPng(file = paste0("fig/hypermutatorBrainRemovedAUC", plotEnding),
          width = 4, height = 4, expr = expression({
            par(mar = c(4,8,2,0.5))
            AUCcollection = c("Original"= AUCorig@y.values[[1]],
                              "Removed from\ntest data" = perfPatients[["brain"]][["TCGA-DU-6392"]]$auc,
                              "Removed from\ntrain data" = AUCfull@y.values[[1]],
                              "Removed from test\nand train data" = AUC@y.values[[1]])
            barplot(rev(AUCcollection-0.5), horiz = T, las = 1, col = 4:1, xlab = "AUC", 
                    xaxt = "n", main = "TCGA-DU-6392", mgp = c(2,0.6,0))
            axis(1, at = axTicks(1), labels = c(0,axTicks(1)[-1]+0.5), mgp = c(2,0.8,0))
            plotrix::axis.break(axis = 1, breakpos = 0.01, style = "gap")
            abline(v=0)
          }))


#####


