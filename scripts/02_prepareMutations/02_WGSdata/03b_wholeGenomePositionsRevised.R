#!/usr/bin/env Rscript
library(GenomicRanges)
library(dplyr)
library(data.table)

# Debug: Print all arguments received
cat("\nRaw arguments received:\n")
print(commandArgs(trailingOnly = TRUE))
cat("\n")

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

# Find argument positions
chrom_idx <- which(args == "--chrom") + 1
bed_idx <- which(args == "--bed") + 1
cmd_idx <- which(args == "--cmd") + 1

# Validate we have all arguments
if(length(args) < 6 || any(c(chrom_idx, bed_idx, cmd_idx) > length(args))) {
  stop("Usage: Rscript script.R --chrom <chrom> --bed <bed_file> --cmd <cmd_file>")
}

# Extract values
chrom <- args[chrom_idx]
bed_file <- args[bed_idx]
cmd_file <- args[cmd_idx]

cat("Parsed arguments:\n")
cat("Chromosome:", chrom, "\n")
cat("BED file:", bed_file, "\n")
cat("CMD file:", cmd_file, "\n\n")
if(!chrom %in% paste0("chr", 1:22)) stop("Chromosome must be chr1-chr22")

  Muts <- fread(bed_file, sep = '\t', header = FALSE, col.names = c('chr','pos','ref','alt'))
  context <- fread(cmd_file, sep = '\t', header = FALSE, col.names = c('context'))


# Read input files
#Muts <- fread(bed_file, sep = '\t', header = FALSE,col.names = c('chr','pos','ref','alt'))
#context <- fread(cmd_file, sep = '\t', header = FALSE,col.names = c('context'))

# Process data
Muts$context <- context$context
Muts$ref <- substr(Muts$context, 3, 3)
Muts$mutated = NA
Muts <- as.data.frame(Muts)

Muts <- Muts %>% 
  filter(ref %in% c("A", "T", "G", "C"))

output_dir <- "data/MutTables/WholeGenomeData/"

# RData output
save(Muts, file = file.path(output_dir, paste0("Muts_WithContext_", chrom, ".RData")))


#ids = do.call(paste,c(Muts2, sep = "_"))
#ids = apply(Muts, 1, paste, collapse="_")
#MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
       #         ids)


#Create IDs with all fields (including NAs) but without extra spaces
ids <- paste0(Muts$chr, "_", Muts$pos, "_", Muts$ref, "_", Muts$alt, "_", Muts$context, "_", Muts$mutated)

# Create BED format (0-based start, 1-based end)
MutsBed <- data.frame(
  chr = Muts$chr,
  start = as.integer(Muts$pos - 1),
  end = as.integer(Muts$pos),
  id = ids,
  stringsAsFactors = FALSE
)



# Save chromosome-specific output
fwrite(MutsBed, file = file.path(output_dir, paste0("MutsResult_", chrom, ".bed")), 
       sep = '\t', col.names = FALSE, quote = FALSE)

print(paste("Saved clean BED for", unique(MutsBed$chr)))
