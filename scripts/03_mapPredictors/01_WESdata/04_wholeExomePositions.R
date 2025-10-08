args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
tissue = tissues[args]
print(tissue)
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(readxl)
source("lib/dataMapping.R")
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)

print("mapping data")
tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
                sheet="allTissues", col_names=T)
tab$NA. = NULL
tab[tab == "NA"] = NA
tab = tab[,c(colnames(tab)[1:9],tissue)]
tab = tab[!is.na(tab[,tissue]),]
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


dumpVar = sapply(1:50, function(i){
  print(paste0("processing part ", i))
  if(file.exists(paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,"_", tissue, "_mapped.RData"))){
    return(NA)
  }
  pred = mapPredictors(x=tab,
                       posFile=paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,".bed"))
  load(paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,".RData"))
  data = list(meta = tab, pred = pred, muts = subMuts)
  save(data, file = paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,"_", tissue, "_mapped.RData"))
  rm(pred, subMuts, data); gc()
  return(NA)
})
  
