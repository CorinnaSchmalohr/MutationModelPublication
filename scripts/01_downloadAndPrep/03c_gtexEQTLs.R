# tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
#             "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)


# gtex  p-value#####
tissue2tissueName = c("brain"="brain_combined",
                      "breast"= "Breast_Mammary_Tissue", #
                      "colon" = "Colon_Sigmoid",
                      "esophagus" = "Esophagus_Mucosa",
                      "kidney" = "Kidney_Cortex", 
                      "liver" = "Liver", 
                      "lung" = "Lung",
                      "ovary" = "Ovary", 
                      "prostate" = "Prostate", 
                      "skin" = "Skin_Sun_Exposed_Lower_leg")
# write each tissue file as a bed file 
dumpVar = sapply(tissue2tissueName, function(tissue){
  print(tissue)
  load(paste0("data/predictors/GTEx_eQTL/", tissue,".RData")) # gr
  colnames(mcols(gr)) = "score" # otherwise the score column is not exported
  gr = gr[width(gr) == 1]
  export(gr, con = paste0("temp/GTEx_eQTL_",tissue,".bedgraph"))
  return(NULL)
})
# combine all  elements and read in R
unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                    paste(" temp/GTEx_eQTL_", tissue2tissueName, ".bedgraph", 
                          sep = "", collapse = " "))
comb = read.table(text = system(unionbedg, intern = T))  
file.remove(paste0("temp/GTEx_eQTL_", tissue2tissueName, ".bedgraph"))
# convert into a granges object
val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
val[rowSums(is.na(comb[,-(1:3)])) == (ncol(comb)-3)] = NA
gr = GRanges(seqnames=comb$V1,
             ranges=IRanges(start=comb$V2+1, end = comb$V3), 
             val = val)
gr = gr[!is.na(gr$val)]
save(gr, file = paste0("data/predictors/GTEx_eQTL/allTissues.RData"))
#####



# GTEx slope #####
tissue2tissueName = c("brain"="brain_slopes_combined",
                      "breast"= "Breast_Mammary_Tissue_slope", #
                      "colon" = "Colon_Sigmoid_slope",
                      "esophagus" = "Esophagus_Mucosa_slope",
                      "kidney" = "Kidney_Cortex_slope", 
                      "liver" = "Liver", 
                      "lung" = "Lung_slope",
                      "ovary" = "Ovary_slope", 
                      "prostate" = "Prostate_slope", 
                      "skin" = "Skin_Sun_Exposed_Lower_leg_slope")
# write each tissue file as a bed file 
dumpVar = sapply(tissue2tissueName, function(tissue){
  print(tissue)
  load(paste0("data/predictors/GTEx_eQTL/", tissue,".RData")) # gr
  gr = sort(gr)
  gr = gr[width(gr) == 1]
  colnames(mcols(gr)) = "score" # otherwise the score column is not exported
  export(gr, con = paste0("temp/GTEx_eQTL_",tissue,".bedgraph"))
  return(NULL)
})
# combine all GRL elements and red in R
unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                    paste(" temp/GTEx_eQTL_", tissue2tissueName, ".bedgraph", 
                          sep = "", collapse = " "))
comb = read.table(text = system(unionbedg, intern = T))  
file.remove(paste0("temp/GTEx_eQTL_", tissue2tissueName, ".bedgraph"))
# convert into a granges object
val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
val[rowSums(is.na(comb[,-(1:3)])) == (ncol(comb)-3)] = NA
gr = GRanges(seqnames=comb$V1,
             ranges=IRanges(start=comb$V2+1, end = comb$V3), 
             val = val)
gr = gr[!is.na(gr$val)]
save(gr, file = paste0("data/predictors/GTEx_eQTL/allTissues_slope.RData"))
#####

