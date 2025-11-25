tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)


# gtex  #####
tissue2tissueName = c("brain"="brain_combined",
                      "breast"= "BreastMammaryTissue", #
                      "colon" = "ColonSigmoid",
                      "esophagus" = "EsophagusMucosa",
                      "kidney" = "KidneyCortex", 
                      "liver" = "Liver", 
                      "lung" = "Lung",
                      "ovary" = "Ovary", 
                      "prostate" = "Prostate", 
                      "skin" = "SkinSunExposedLowerleg")

load(paste0("data/predictors/GTEx_expression/",
            tissue2tissueName[1], ".RData"))
expr = sapply(tissue2tissueName, function(tissue){
  load(paste0("data/predictors/GTEx_expression/",
              tissue, ".RData"))
  gr$val
})
gr$val = rowMeans(expr, na.rm = T)
save(gr, file = "data/predictors/GTEx_expression/allTissues.RData")
#####


# cancer expression #####
tissue2tissueName = c("brain"="lgg",
                      "breast"= "brca", #
                      "colon" = "coadread",
                      "esophagus" = "esca",
                      "kidney" = "kirc", 
                      "liver" = "lihc", 
                      "lung" = "luad",
                      "ovary" = "ov", 
                      "prostate" = "prad", 
                      "skin" = "skcm")
load(paste0("data/predictors/cancerExpression/",
            tissue2tissueName[1], ".RData"))
expr = sapply(tissue2tissueName, function(tissue){
  load(paste0("data/predictors/cancerExpression/",
              tissue, ".RData"))
  gr$val
})
gr$val = rowMeans(expr, na.rm = T)
save(gr, file = "data/predictors/cancerExpression/allTissues.RData")
#####




# pancan healthy
tissue2tissueName = c("breast"= "brca", #"brain"="lgg",
                      "colon" = "coadread",
                      "esophagus" = "esca",
                      "kidney" = "kirc", 
                      "liver" = "lihc", 
                      "lung" = "luad",
                      # "ovary" = "ov", 
                      "prostate" = "prad" #, "skin" = "skcm"
                      )
load(paste0("data/predictors/TCGAnormals_expression/",
            tissue2tissueName[1], ".RData"))
expr = sapply(tissue2tissueName, function(tissue){
  load(paste0("data/predictors/TCGAnormals_expression/",
              tissue, ".RData"))
  gr$val
})
gr$val = rowMeans(expr, na.rm = T)
gr$val[is.nan(gr$val)] = NA
save(gr, file = "data/predictors/TCGAnormals_expression/allTissues.RData")
