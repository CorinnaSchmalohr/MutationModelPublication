args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
muttissue = tissues[args]
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
source("lib/dataMapping.R")
library(readxl)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
print(muttissue)
dir.create("data/MutTables/exomeTrainDataCrossTissue", showWarnings = F)

print("mapping cross-tissue data")
# cross-tissue: for each tissue, prepare predictors from all other tissues

for(predtissue in tissues[tissues != muttissue]){
  print(predtissue)
  tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
                  sheet="allTissues", col_names=T)
  tab$NA. = NULL
  tab[tab == "NA"] = NA
  tab = tab[,c(colnames(tab)[1:9],predtissue)]
  tab = tab[!is.na(tab[,predtissue]),]
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
  
  pred = mapPredictors(x=tab, 
                       posFile=paste0("data/MutTables/exomeTrainData/", 
                                      muttissue, "_Muts.bed"))
  load(paste0("data/MutTables/exomeTrainData/", muttissue, "_Muts.RData"))
  data = list(meta = tab, pred = pred, muts = Muts)
  save(data, file = paste0("data/MutTables/exomeTrainDataCrossTissue/",
                           muttissue, "With",predtissue,"Preds", ".RData"))
  dat = cbind(data$pred, mutated = as.factor(data$muts$mutated))
  datchroms = data$muts$chr
  dat = as.data.frame(dat)
  save(dat,datchroms, file=paste0("data/MutTables/exomeTrainDataCrossTissue/", 
                                  muttissue, "With",predtissue,
                                  "Preds_processed", ".RData"))
}
cat("\n")



