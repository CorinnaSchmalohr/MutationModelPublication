tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)

# get tissues for which we have data
tissues = tissues[file.exists(paste0("data/predictors/DNAseENCODE/", 
                                     tissues,".RData"))]
# write each tissue file as a bed file 
dumpVar = sapply(tissues, function(tissue){
  load(paste0("data/predictors/DNAseENCODE/", tissue,".RData")) # gr
  seqlevels(gr) = paste0("chr", c(1:22, "X", "Y"))
  gr = sort(gr)
  colnames(mcols(gr)) = "score" # otherwise the score column is not exported
  export(gr, con = paste0("temp/DNAse_",tissue,".bed.gz"))
  return(NULL)
})
# merge overlapping intervals
fixOverlapping = paste0("bedtools merge  -c 5 -o mean -d -1 ",
                        "-i temp/DNAse_", tissues, ".bed.gz ",
                        "> temp/DNAse_", tissues, "_merged.bedgraph ")
cmd=paste(fixOverlapping, collapse = "\n")
system(cmd)
# combine all GRL elements and red in R
unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                    paste(" temp/DNAse_", tissues, "_merged.bedgraph", 
                          sep = "", collapse = " "))
comb = read.table(text = system(unionbedg, intern = T))  
file.remove(c(paste0("temp/DNAse_", tissues, ".bed.gz"), 
              paste0("temp/DNAse_", tissues, "_merged.bedgraph")))
# convert into a granges object
val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
gr = GRanges(seqnames=comb$V1,
             ranges=IRanges(start=comb$V2+1, end = comb$V3), 
             val = val)
save(gr, file = paste0("data/predictors/DNAseENCODE/allTissues.RData"))
