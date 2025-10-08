library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
source("scripts/05_analysis/00_NamesAndColors.R")
rarePredictors = c("ZBTB33_100bp", "YY1_100bp", "TAF1_100bp", "SP1_100bp",
                   "RXRA_100bp", "REST_100bp", "RAD21_100bp",
                   "NR2F2_100bp", "MAX_100bp", "JUND_100bp", "HNF4G_100bp",
                   "HNF4A_100bp", "GABPA_100bp", "FOXA2_100bp", 
                   "FOXA1_100bp", "EGR1_100bp", "ATF3_100bp")
predictorOrder = read.table("scripts/04_trainmodels/predictorOrder_glm.txt")[,1]


# Generate models #####
dumpVar =   lapply(tissues, function(tissue){
  print(tissue)
    
  # prepare data for this tissue ######
  full_dat = sapply(tissues[tissues != tissue], function(t){
    tab = read_xlsx("data/rawdata/dataMapping.xlsx", sheet = t, col_names =T)
    features = tab[tab$sourceTissue == t,]$abbreviation
    load(paste0("data/procData/traindata_crosstissue/traindata_crosstissue_",tissue,"With",t,"Preds_processed.RData"))
    dat = dat[features]
    return(dat)
  })
  full_dat = do.call(cbind, full_dat)
  
  
  tab = read_xlsx("data/rawdata/dataMapping.xlsx", sheet = tissue, col_names =T)
  features = tab[tab$sourceTissue == tissue,]$abbreviation
  load(paste0("data/procData/traindata/traindata_processed_", tissue, ".RData"))
  colnames(dat[features]) = paste0(tissue,".",colnames(dat[features]))
  setnames(dat, old = colnames(dat[features]), new =  paste0(tissue,".",colnames(dat[features]))) # rename only tissue specific features
  
  full_dat = cbind(dat, full_dat)
  chroms = unique(datchroms)
  rm(dat, tab)
  #####
  
  # glm #####
  print("training glm")
  logR = glm(formula = mutated ~ ., data = full_dat, family = binomial(link = "logit"))
  save(logR, file = paste0("data/rdata/GLMmodel/", tissue, "_allFeatures_selection.RData"))
  
  print("p-values")
  pvals = coef(summary(logR))[,4][-1]
  sigFreatures = names(pvals[pvals < 0.05])
  
  # See how many are picked from the source tissue
  list_sigF = strsplit(x = sigFreatures, split = ".", fixed = TRUE)
  
  selection = table(unlist(list_sigF))
  selection_tissues = selection[names(selection) %in% tissues]
  missing_tissues = tissues[!tissues %in% names(selection_tissues)]
  selection_tissues[missing_tissues] = 0
  selection_tissues = selection_tissues[tissues]
  
  
  
  # Plots
  c = rep("grey50", length(selection_tissues))
  names(c) = names(selection_tissues)
  c[tissue] = tissueCols[tissue]
  png(paste0("fig/modelEvaluation/FeatureSource_SelfSelection_",tissue,".png"),
      width=1500, height=800, res=150)
  plt = barplot(as.vector(selection_tissues), xlab = "Source tissue of selected features", ylab = "Number of selections", 
                names=t2T[names(selection_tissues)], las = 1, col = c)
  text(x = plt, y = selection_tissues-0.5, label = selection_tissues)  
  dev.off()
  #####
  
  print("done with this tissue")
})
#####

# Plot overview of feature sources for all tissues #####
feature_countT = setNames(data.frame(matrix(ncol = length(tissues), nrow = 0)), tissues)
for(tissue in tissues){
  load(file = paste0("data/rdata/GLMmodel/", tissue, "_allFeatures_selection.RData"))
  
  print(tissue)
  pvals = coef(summary(logR))[,4][-1]
  sigFreatures = names(pvals[pvals < 0.05])
  
  # See how many are picked from the source tissue
  list_sigF = strsplit(x = sigFreatures, split = ".", fixed = TRUE)
  for(l in list_sigF){
    if(length(l) > 1){
      t = l[1]
      f = l[2]
      if(f %in% rownames(feature_countT)){
        feature_countT[f, t] = feature_countT[f, t] + 1
      } else {
        feature_countT[f, t] = 1
      }
    }
  }
  feature_countT[is.na(feature_countT)] <- 0
}
# Sort data 
feature_countT = feature_countT[na.omit(match(predictorOrder, rownames(feature_countT))),]

feature_countT_noTF = feature_countT[!rownames(feature_countT) %in% rarePredictors,]
rownames(feature_countT_noTF) = p2P[rownames(feature_countT_noTF)]
rownames(feature_countT) = p2P[rownames(feature_countT)]

feature_countT_percentage = apply(feature_countT, 1, function(x){x*100/sum(x,na.rm=T)})
feature_countT_percentage_noTF = apply(feature_countT_noTF, 1, function(x){x*100/sum(x,na.rm=T)})

png(paste0("fig/modelEvaluation/FeatureSource_SelfSelection.png"),
    width=1500, height=1000, res=150)
par(mfrow = c(3,1), mar = c(1,3,0,0))
barplot(as.matrix(t(feature_countT)), border="white", xlab="", col = tissueCols,  xaxt = "n", las = 1)
title(ylab="Count", line=2, cex.lab=1.2)
legend("topright", rownames(feature_countT_percentage), fill = tissueCols[rownames(feature_countT_percentage)], cex = 1.2)
barplot(feature_countT_percentage, border="white", xlab="Feature", col = tissueCols, las = 2)
title(ylab="Percent %", line=2, cex.lab=1.2)
par(mar = c(2,3,10,0))
barplot(t(rowSums(feature_countT_percentage)/ncol(feature_countT_percentage)), border="white", names.arg = t2T)
title(ylab="Percent %", line=2, cex.lab=1.2)
dev.off()

# Ignore rare predictors 
png(paste0("fig/modelEvaluation/FeatureSource_SelfSelection_noTF.png"),
    width=1500, height=1000, res=150)
par(mfrow = c(3,1), mar = c(1,3,0,0))
barplot(as.matrix(t(feature_countT_noTF)), border="white", xlab="", col = tissueCols,  xaxt = "n", las = 1)
title(ylab="Count", line=2, cex.lab=1.2)
legend("topleft", rownames(feature_countT_percentage_noTF), fill = tissueCols[rownames(feature_countT_percentage_noTF)],horiz = T, bty = "n")
barplot(feature_countT_percentage_noTF, border="white", xlab="Feature", col = tissueCols, las = 2)
title(ylab="Percent %", line=2, cex.lab=1.2)
par(mar = c(2,3,10,0))
barplot(t(rowSums(feature_countT_percentage_noTF)/ncol(feature_countT_percentage_noTF)), border="white", names.arg = t2T)
title(ylab="Percent %", line=2, cex.lab=1.2)
dev.off()
#####



# Rank the predictors by p-value and check, which tissue sources are in the top ranked positions #####
feature_rankingT = setNames(data.frame(matrix(ncol = length(tissues), nrow = 0)), tissues)
first_rank = c()
for(tissue in tissues){
  print(tissue)
  load(file = paste0("data/rdata/GLMmodel/", tissue, "_allFeatures_selection.RData"))
  ranked_pvals = sort(coef(summary(logR))[,4][-1])
  sigFreatures = names(ranked_pvals[ranked_pvals < 0.05])
  sigFreatures = sigFreatures[grep("\\.", sigFreatures)] # select only the tissue-specific predictors 

  ranked_list = strsplit(x = sigFreatures, split = "\\.")
  ranked_list = as.data.frame(do.call(rbind, ranked_list))
  colnames(ranked_list) = c("tissue", "feature")
  ranked_list["rank"] = 1:nrow(ranked_list)
  
  # List the tissue source of the highest ranked feature for each model 
  first_rank = c(first_rank, ranked_list$tissue[1])
  
  # Check the highest ranked source of features that were selected multiple times coming from different tissues
  multiple_features = names(which(table(ranked_list$feature) > 1))
  #ranked_list[ranked_list$feature %in% multiple_features,]
  for(f in multiple_features){
    feature_rankingT[f, tissue] = ranked_list[ranked_list$feature == f,]$tissue[1] # Source of highest ranked feature
  }
}

names(first_rank) = tissues
same_sources = names(first_rank) == first_rank
png(paste0("fig/modelEvaluation/FeatureSource_SelfSelection_highestRank.png"),
    width=1000, height=500, res=130)
barplot(c(sum(same_sources), sum(!same_sources)), xlab = "Source of the highest ranked feature", las = 1, 
        names.arg = c("Feature source = Mutation source", "Feature source != Mutation source"	))
text(x = 0.7, y = sum(same_sources)-1, label = paste(first_rank[same_sources],collapse=' '))  
text(x = 1.9, y = sum(!same_sources)-1, label = paste(names(first_rank[!same_sources]),collapse=' '))  
text(x = 1.9, y = sum(!same_sources)-1.5, label = paste(first_rank[!same_sources],collapse=' '))  
dev.off()


same_source_multF_percentage = sapply(colnames(feature_rankingT), function(t){
  # Counts, how often the source of the highest ranked feature was the "correct" one, if the feature was sig. multiple time in the model
  count_sameSource = sum(feature_rankingT[t] == t, na.rm = TRUE)
  return(count_sameSource*100/sum(!is.na(feature_rankingT[t])))
})
barplot(same_source_multF_percentage, ylab = "Percent of 'correct' source")
