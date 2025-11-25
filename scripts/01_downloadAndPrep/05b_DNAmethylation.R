

tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)
library(parallel)

# get tissues for which we have data
tissues = tissues[file.exists(paste0("data/predictors/methylationGSE186458/", 
                                     tissues, ".RData"))]



GRL = sapply(paste0("chr", 1:22), function(chrom){
  cat(chrom, ' ')
  GRLchrom  = lapply(tissues, function(tissue){
    load(paste0("data/predictors/methylationGSE186458/", 
                tissue, ".RData")) # gr
    gr = gr[seqnames(gr) == chrom]
    return(gr)
  }) 
  names(GRLchrom) = tissues
  gr <- do.call(c, c(GRLchrom, use.names = FALSE))
  mcols(gr) <- NULL
  gr <- unique(sort(gr))
  # get dataframe of signal value for each dataset in the output GRanges
  idx <- mclapply(GRLchrom, function(x) which(gr %in% x), mc.cores = 4)
  counts <- mcmapply(function(dat, idx) {
    out <- rep.int(NA, length(gr))
    out[idx] <- mcols(dat)[["score"]]
    out
  }, GRLchrom, idx, mc.cores = 8, SIMPLIFY = TRUE)
  gr$score = rowSums(counts, na.rm = T) / ncol(counts)
  return(gr)
})
gr = do.call(c, c(GRL, use.names = FALSE))
save(gr, file = paste0("data/predictors/methylationGSE186458/allTissues.RData"))
