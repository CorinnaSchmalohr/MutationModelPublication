library(GenomicRanges)
library(rtracklayer)
tissues =c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
           "ovary", "prostate", "skin")

# interaction frequencies #####
GRLints = sapply(tissues, function(tissue){
  print(tissue)
  load(paste0("data/predictors/HiCENCODE/", tissue, "_ints.RData")) # GRLints
  ids = names(GRLints)
  # write each element as a bed file 
  dumpVar = sapply(ids, function(id){
    GRLelement = GRLints[[id]]
    seqlevels(GRLelement) = paste0("chr", c(1:22, "X", "Y"))
    GRLelement = sort(GRLelement)
    colnames(mcols(GRLelement)) = "score" # otherwise the score column is not exported
    export(GRLelement, con = paste0("temp/",id,".bed.gz"))
    return(NULL)
  })
  # merge overlapping intervals
  fixOverlapping = paste0("bedtools merge  -c 5 -o mean -d -1 ",
                          "-i temp/", ids, ".bed.gz ",
                          "> temp/", ids, "_merged.bedgraph ")
  cmd=paste(fixOverlapping, collapse = "\n")
  system(cmd)
  # combine all GRL elements and red in R
  unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                      paste(" temp/", ids, "_merged.bedgraph", sep = "", collapse = " "))
  comb = read.table(text = system(unionbedg, intern = T))  
  file.remove(c(paste0("temp/", ids, ".bed.gz"), 
                paste0("temp/", ids, "_merged.bedgraph")))
  # convert into a granges object
  comb = comb[comb$V1 %in% paste0("chr", 1:22),]
  val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
  gr = GRanges(seqnames=comb$V1,
               ranges=IRanges(start=comb$V2+1, end = comb$V3), 
               score = val)
  return(gr)
}, simplify = F)
GRLints = GRangesList(GRLints)
save(GRLints, file = paste0("data/predictors/HiCENCODE/allTissues_ints.RData"))
#####

# compartments #####
GRLcomp = sapply(tissues, function(tissue){
  print(tissue)
  if(!file.exists(paste0("data/predictors/HiCENCODE/", tissue, "_comp.RData"))){
    return(NULL)
  }
  load(paste0("data/predictors/HiCENCODE/", tissue, "_comp.RData")) # GRLcomp
  ids = names(GRLcomp)
  if(length(ids) < 2){
    return(GRLcomp[[1]])
  }
  # write each element as a bed file 
  dumpVar = sapply(ids, function(id){
    GRLelement = GRLcomp[[id]]
    seqlevels(GRLelement) = paste0("chr", c(1:22, "X", "Y"))
    GRLelement = sort(GRLelement)
    colnames(mcols(GRLelement)) = "score" # otherwise the score column is not exported
    export(GRLelement, con = paste0("temp/",id,".bed.gz"))
    return(NULL)
  })
  # merge overlapping intervals
  fixOverlapping = paste0("bedtools merge  -c 5 -o mean -d -1 ",
                          "-i temp/", ids, ".bed.gz ",
                          "> temp/", ids, "_merged.bedgraph ")
  cmd=paste(fixOverlapping, collapse = "\n")
  system(cmd)
  # combine all GRL elements and red in R
  unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                      paste(" temp/", ids, "_merged.bedgraph", sep = "", collapse = " "))
  comb = read.table(text = system(unionbedg, intern = T))  
  file.remove(c(paste0("temp/", ids, ".bed.gz"), 
                paste0("temp/", ids, "_merged.bedgraph")))
  # convert into a granges object
  comb = comb[comb$V1 %in% paste0("chr", 1:22),]
  val = rowSums(comb[,-(1:3)], na.rm = T) / (ncol(comb)-3)
  gr = GRanges(seqnames=comb$V1,
               ranges=IRanges(start=comb$V2+1, end = comb$V3), 
               score = val)
  return(gr)
}, simplify = F)
GRLcomp = GRLcomp[!sapply(GRLcomp, is.null)]
GRLcomp = GRangesList(GRLcomp)
save(GRLcomp, file = paste0("data/predictors/HiCENCODE/allTissues_comp.RData"))
#####
