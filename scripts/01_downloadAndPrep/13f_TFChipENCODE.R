tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
library(GenomicRanges)
library(rtracklayer)

# get tissues and Chipseq target combos for which we have data
dataAvail = sapply(tissues, function(tissue){
  meta = read.table(paste0("data/rawdata/TFChipENCODE/", tissue, "/metadata.tsv"),
                     header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type == "bed", ]
  return(unique(meta$Experiment.target))
})
targets = unique(unlist(dataAvail))

dumpVar = sapply(targets, function(target){
  print(target)
  # get available tissues
  tissueAvail = file.exists(paste0("data/predictors/TFChipENCODE/",
                                   target, "_", tissues,".RData"))
  tempTissues = tissues[tissueAvail]
  # write each tissue file as a bed file 
  dumpVar2 = sapply(tempTissues, function(tissue){
    load(paste0("data/predictors/TFChipENCODE/",
                target, "_", tissue,".RData")) # gr
    seqlevels(gr) = paste0("chr", c(1:22, "X", "Y"))
    gr = sort(gr)
    colnames(mcols(gr)) = "score" # otherwise the score column is not exported
    export(gr, con = paste0("temp/TFChip_", target, "_", tissue, ".bed.gz"))
    return(NULL)
  })
  # merge overlapping intervals
  fixOverlapping = paste0("bedtools merge  -c 5 -o mean -d -1 ",
                          "-i temp/TFChip_", target, "_", tempTissues, ".bed.gz ",
                          "> temp/TFChip_", target, "_", tempTissues, "_merged.bedgraph ")
  cmd=paste(fixOverlapping, collapse = "\n")
  system(cmd)
  # combine all GRL elements and read in R
  unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                      paste(" temp/TFChip_", target, "_", tempTissues, "_merged.bedgraph", 
                            sep = "", collapse = " "))
  comb = read.table(text = system(unionbedg, intern = T))  
  file.remove(c(paste0("temp/TFChip_", target, "_", tempTissues, ".bed.gz"), 
                paste0("temp/TFChip_", target, "_", tempTissues, "_merged.bedgraph")))
  # convert into a granges object
  val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
  gr = GRanges(seqnames=comb$V1,
               ranges=IRanges(start=comb$V2+1, end = comb$V3), 
               val = val)
  save(gr, file = paste0("data/predictors/TFChipENCODE/", target, "_", "allTissues.RData"))
})




