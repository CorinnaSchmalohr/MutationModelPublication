args = as.numeric(commandArgs(trailingOnly=T))
tissues = c("brain","breast", "colon","esophagus","kidney","liver", "lung","ovary",
            "prostate", "skin")
tissue = tissues[args]
print(tissue)
library(ranger)
nThreads = 24
nPermutations = 10000
maxData = 50000

# prepare data #####
#load(paste0(paste0("data/procData/traindata/traindata_processed_",
#                   tissue, ".RData")))
#chroms = unique(datchroms)
#if(nrow(dat) > maxData){
#   samp = sample(1:nrow(dat), maxData, replace=F)
#   dat = dat[samp,]
#   save(samp, dat,file = paste0("data/rdata/RFmodel/",tissue, "_subDataForPvals.RData"))
#}
# Use downsmapled data from following scripts 
load(paste0("data/procData/traindata/traindata_processed_",tissue, "_subsampled.RData"))
######


# grow forest with impurity_corrected ####
#print("growing forests with impurity_corrected")
#temp = lapply(chroms, function(cr){
#   cat(cr, ' ')
#   trainData = dat[datchroms != cr,]
#   rf = ranger(mutated ~ ., data = trainData, importance = 'impurity_corrected',
#               write.forest = T, seed = 1234, num.threads =  nThreads,
#               respect.unordered.factors = 'partition',
#               probability = T, verbose=T)
#   save(rf, file = paste0("data/rdata/RFmodel/", tissue, "_", 
#                          cr, "_impcorr.RData"))
#   return(NULL)
#})
# and extract these importances
#print("extracting importances")
#imp = sapply(chroms, function(cr){
#   cat(cr, ' ')
   #load(paste0("data/rdata/RFmodel/", tissue, "_", cr, "_impcorr.RData"))
#   load(paste0("data/rdata/RFmodel/", tissue, "_subsampled_", cr, ".RData"))
#   return(rf$variable.importance)
#})
#save(imp, file = paste0("data/rdata/RFmodel/", tissue,
#                        "_importances_impcorr.RData"))
######


# calculate p-values based on whole data with impurity_corrected #####
print("compute rf on whole data")
rf = ranger(mutated ~ ., data = dat, importance='impurity_corrected',
            write.forest = T, seed = 1234, num.threads = nThreads,
            respect.unordered.factors = 'partition',
            probability = T, verbose=T)
save(rf, file = paste0("data/rdata/RFmodel/", tissue,
                       "_wholeDataRF.RData"))
print("permutations")
pvals = importance_pvalues(rf, method="altmann", num.permutations=nPermutations,
                           mutated ~ ., data = dat,
                           seed = 1234, num.threads = nThreads,
                           respect.unordered.factors = 'partition',
                           scale.permutation.importance = T,
                           probability = T, verbose=T)
save(pvals, file = paste0("data/rdata/RFmodel/", tissue,
                          "_pvals.RData"))
print("done")
#####
