library(ranger)
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
chr_num <- as.integer(args[1])
part_num <- as.integer(args[2])  # New argument for part number
tissues <- c("esophagus")
#tissues <- c("brain","breast","esophagus", "kidney", "liver", "ovary", "prostate", "skin")

nThreads <- 56

for (tissue in tissues) {
  cat("\n==== Processing chromosome", chr_num, "part", part_num, "for tissue:", tissue, "====\n")
  
  # Format part number with leading zeros
  part_num_padded <- sprintf("%03d", part_num)
  
  input_file <- paste0("data/MutTables/WholeGenomeData/partial/esophagus/",
                       tissue, "_MutsResult_chr", chr_num, "_part", part_num_padded, ".RData")
  
  output_file <- paste0("data/Modeling/WholeGenomeData/RFnew/esophagus/",
                        tissue, "_MutsResult_chr", chr_num, "_part", part_num_padded, "_RFpredictions.RData")
  
  model_file <- paste0("data/Modeling/WholeGenomeData/RF/",
                       tissue, "_finalModel.RData")
  
  # Skip conditions
  if (!file.exists(input_file)) {
    cat("SKIP: Input file not found:", input_file, "\n")
    next
  }
  if (file.exists(output_file)) {
    cat("SKIP: Output already exists:", output_file, "\n")
    next
  }
  if (!file.exists(model_file)) {
    cat("SKIP: Model file not found:", model_file, "\n")
    next
  }
  
  # Load data
  cat("Loading data from:", input_file, "\n")
  load(input_file)
  
  # Load model
  cat("Loading model from:", model_file, "\n")
  load(model_file)   
  
  # Make predictions
  cat("Making predictions using model features...\n")
  
  testDat <- partial_data$pred
  yhat <- predict(rf, testDat, num.threads = nThreads)
  
  predictions <- data.frame(partial_data$muts, prediction = yhat$predictions[,2])
  
  save(predictions, file = output_file)
  
  cat("SUCCESS: Processed chromosome", chr_num, "part", part_num, "for tissue:", tissue, "\n")
}