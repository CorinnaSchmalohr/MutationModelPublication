library(GenomicRanges)
# load exonic regions
load("data/processedData/codExons_noUTR_filtered.RData")



dir.create("data/MutTables/WholeExomeData", showWarnings = F)
# get all positions in the exome
Muts = lapply(paste0("chr", 1:22), function(cr){
  print(cr)
  chrexons = codExons_filtered[codExons_filtered$chr == cr,]
  pos = do.call(c,apply(chrexons,1, function(x){
    x["start"]:x["end"]
  }))
  pos = unique(pos)
  m = data.frame(chr = cr, pos = pos, ref = NA, alt = NA)
  return(m)
})
Muts = do.call(rbind, Muts)

# get context 
bed = data.frame(Muts$chr, as.integer(Muts$pos-3), as.integer(Muts$pos+2))
write.table(bed, file="temp/TPs_bed.bed", quote=F, col.names=F, row.names=F, sep="\t")
cmd = paste0("bedtools getfasta -fi data/rawdata/genome/GRCh37.primary_assembly.genome.fa ",
             "-bed temp/TPs_bed.bed ",
             "-name -tab | cut -f 2")   
context = system(cmd, intern = T)
Muts$context = context
Muts$ref = substr(Muts$context, 3,3)
Muts$mutated = NA

# save 
save(Muts, file = "data/MutTables/WholeExomeData/exomeMuts.RData")
ids = do.call(paste,c(Muts, sep = "_"))
MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                ids)
write.table(MutsBed, file = "data/MutTables/WholeExomeData/exomeMuts.bed",
            col.names = F, row.names = F, sep = "\t", quote = F)
rm(ids, MutsBed, bed)

# save in parts because of memory problems otherwise
splitVar = gl(n=50, k= round(nrow(Muts)/50), length = nrow(Muts))
splitMuts = split(Muts, splitVar)
dumpVar = sapply(1:50, function(i){
  Muts = splitMuts[[i]]
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                  ids)
  save(Muts, file = paste0("data/MutTables/WholeExomeData/exomeMuts_part", i,".RData"))
  write.table(MutsBed, file = paste0("data/MutTables/WholeExomeData/exomeMuts_part", i, ".bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)
})
#####



# create binned version of exome #####
codExonsGR = GRanges(codExons_filtered)
codExonsGR_merged <- reduce(codExonsGR)
bins100 = slidingWindows(codExonsGR_merged, width = 100, step = 100)
bins100 = unlist(bins100)
bins1kb = slidingWindows(codExonsGR_merged, width = 1000, step = 1000)
bins1kb = unlist(bins1kb)
hist(width(bins100))
hist(width(bins1kb))
save(bins100, bins1kb, file = "data/processedData/exomeBins.RData")
#####