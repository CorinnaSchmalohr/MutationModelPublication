library(GenomicRanges)
library(rtracklayer)
meta = read.table("data/rawdata/replication/metadata.tsv",
                  sep = "\t", header = T)

# Combine all Peaks and Valleys #####
dumpVar = sapply(c("Peaks", "Valleys"), function(view){
  print(view)
  ids = meta$tableName[meta$view == view]
  # sometimes there are overlapping intervals. Since we don't use the score, simply merge them
  # -c and -o options are necessary so that we get a fourth column
  fixOverlapping = paste0("bedtools merge -c 5 -o distinct -d -1 ",
                          "-i data/rawdata/replication/", ids, ".bed.gz",
                          " > temp/", ids, "_merged.bedgraph ")
  cmd=paste(c(fixOverlapping), collapse = "\n")
  system(cmd)
  print("done with overlapping")
  # combine multiple bedgraph files into one.
  unionbedg =  paste0("bedtools unionbedg -filler NA -i ", 
                      paste("temp/", ids, "_merged.bedgraph", sep = "", collapse = " "))
  comb = read.table(text = system(unionbedg, intern = T)) 
  print("done with combining")
  # remove temp files
  file.remove(paste0("temp/", ids, "_merged.bedgraph"))
  # reduce to autosomes
  comb = comb[comb$V1 %in% paste0("chr", 1:22),]
  # compute how many of the files had a peak there
  val = apply(comb[,-(1:3)],1,function(x){sum(!is.na(x))})/(ncol(comb)-3)
  gr = GRanges(seqnames=comb$V1,
               ranges=IRanges(start=comb$V2+1, end = comb$V3), 
               val = val)
  save(gr, file = paste0("data/predictors/replication/allTissues_",view,".Rdata"))
})
#####

# Combine all WaveSignal files #####
print("now doing WaveSignal")
ids = meta$tableName[meta$view == "WaveSignal"]
# I manually checked, all bigWigs use the same 1000bp windows
gr = import(paste0("data/rawdata/replication/", 
                   ids[1], ".bigWig"))
tempScores = sapply(ids, function(id){
  tempGR = import(paste0("data/rawdata/replication/", 
                     id, ".bigWig"))
  tempGR$score
})
vals= rowMeans(tempScores, na.rm = T)
gr$score = (vals-min(vals))/(max(vals)-min(vals))
save(gr, file = paste0("data/predictors/replication/allTissues_WaveSignal.Rdata"))
#####



