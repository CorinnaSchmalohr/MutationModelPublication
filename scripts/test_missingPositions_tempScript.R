i = 1
tissue = "brain"
load(paste0("data/MutTables/WholeGenomeData/Muts_WithContext_chr",i,".RData")) # Muts
# saving some memory
Muts$alt = NULL 
Muts$context = NULL
Muts$mutated = NULL
gc()
load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_MutsResult_chr",
            i,"_parts_1_50_RFpredictions.RData")) # predictions
RFpred1 = predictions
# part 2
load(paste0("data/Modeling/WholeGenomeData/RF/", tissue, "_MutsResult_chr",
            i,"_parts_51_100_RFpredictions.RData")) # predictions


nrow(Muts)
nrow(RFpred1) + nrow(predictions)

missing = !(Muts$pos %in% RFpred1$pos | Muts$pos %in% predictions$pos)
missing = which(missing)
Muts[missing[1:10],]

