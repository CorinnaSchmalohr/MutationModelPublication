args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
tissue = tissues[args]
print(tissue)

nThreads = 28
.libPaths(new = "/data/public/cschmalo/R-4.1.2/")
library(ranger)


dir.create("data/Modeling/WholeExomeData/RF", showWarnings = F)

load(paste0("data/Modeling/exomeTrainData/RF/", tissue, "_finalModel.RData"))
dumpVar = sapply(1:50, function(i){
   print(i)
   # load data
   load(paste0("data/MutTables/WholeExomeData/exomeMuts_part",i,"_", tissue, "_mapped.RData"))
   testDat = data$pred
   yhat = predict(rf, testDat, num.threads = nThreads)
   
   predictions = data.frame(data$muts,
                            prediction = yhat$predictions[,2])
   save(predictions, 
        file = paste0("data/Modeling/WholeExomeData/RF/exomeMuts_part",
                      i, "_",tissue, "_RFpredictions.RData"))
   return(NA)
})
print("done")
