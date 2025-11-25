library(readxl)

# Tissues 
tissues = c("brain","breast", "colon","esophagus", 
            "kidney", "liver", "lung","ovary", 
            "prostate", "skin") #
WGStissues = c("brain", "breast", "esophagus", "kidney", "liver",
               "ovary", "prostate", "skin")
#0aff99#tissueCols = setNames(rainbow(length(tissues)), tissues)
tissueCols = setNames(c("#ff0000","#ff8700","#ffd300","#a1ff0a",
                        "#0aff99","#0aefff","#147df5","#580aff",
                        "#be0aff", "#FF0099"), tissues)

t2T = setNames(c("Brain","Breast", "Colon","Esophagus", 
                 "Kidney", "Liver", "Lung","Ovary", 
                 "Prostate", "Skin"), tissues)

methodCols = setNames(RColorBrewer::brewer.pal(3, "Dark2"), c("RF", "GLM", "SL"))

# Predictors
# get predictor order
tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
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
save(predictorOrder, file = "data/processedData/predictorOrder.RData")
#####


p2P = setNames(predictorOrder$Name, predictorOrder$abbreviation)
predictorGroups = split(predictorOrder$abbreviation, predictorOrder$Group)
p2G = setNames(predictorOrder$Group, predictorOrder$abbreviation)
rm(tab, predictorOrder)

# Chromosomes
chrCols = setNames(rainbow(22), sort(paste0("chr", 1:22)))

# Bases
#basesCol = c("#f18aad","#ea6759","#f88f58","#f3c65f","#8bc28c","#6667ab") #brewer.pal(n = 6, name = "Dark2")
basesCol = c("#14baeb","#030303","#df2b27","#999999","#a2ca61","#ebc7c3") # Same as used by Alexandrov et al.



# Predictors
# get predictor order
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
predictorOrderWGS = tab[,1:3] # Group, Name, abbreviation
save(predictorOrderWGS, file = "data/processedData/predictorOrderWGS.RData")
#####


p2PWGS = setNames(predictorOrderWGS$Name, predictorOrderWGS$abbreviation)
predictorGroups = split(predictorOrderWGS$abbreviation, predictorOrderWGS$Group)
p2GWGS = setNames(predictorOrderWGS$Group, predictorOrderWGS$abbreviation)
rm(tab, predictorOrderWGS)
