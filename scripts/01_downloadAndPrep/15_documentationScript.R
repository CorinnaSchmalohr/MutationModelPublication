tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")

# TFChipENCODE #####
IDS = do.call(rbind, sapply(tissues, function(tissue){
  meta = read.table(paste0("data/rawdata/TFChipENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type %in% c("bigWig", "bed") , ]
  targets = unique(meta$Experiment.target)
  filesByTargetBW = split(meta$File.accession[meta$File.type == "bigWig"],
                          meta$Experiment.target[meta$File.type == "bigWig"])
  filesByTargetBed = split(meta$File.accession[meta$File.type == "bed"],
                           meta$Experiment.target[meta$File.type == "bed"])
  res = sapply(targets, function(target){
    rbind(cbind(tissue, target, "bigWig", paste(filesByTargetBW[[target]], collapse = ", ")), 
          cbind(tissue, target, "bed", paste(filesByTargetBed[[target]], collapse = ", ")))
  }, simplify = F)
  res = do.call(rbind, res)
  return(res)
}, simplify = F))
# View(IDS)
write.table(IDS, file = "temp/TFChipENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)

metafile = do.call(rbind,sapply(tissues, function(tissue){
  temp = read.table(paste0("data/rawdata/TFChipENCODE/", tissue, "/metafile.txt"))
  return(c(tissue, temp))
}, simplify = F))
write.table(metafile, file = "temp/TFChipENCODE_metafile.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = F)
#####

# histone ENCODE #####
IDS = do.call(rbind, sapply(tissues, function(tissue){
  meta = read.table(paste0("data/rawdata/histoneENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type %in% c("bigWig", "bed") , ]
  
  targets = unique(meta$Experiment.target)
  filesByTargetBW = split(meta$File.accession[meta$File.type == "bigWig"],
                          meta$Experiment.target[meta$File.type == "bigWig"])
  filesByTargetBed = split(meta$File.accession[meta$File.type == "bed"],
                           meta$Experiment.target[meta$File.type == "bed"])
  res = sapply(targets, function(target){
    rbind(cbind(tissue, target, "bigWig", paste(filesByTargetBW[[target]], collapse = ", ")), 
          cbind(tissue, target, "bed", paste(filesByTargetBed[[target]], collapse = ", ")))
  }, simplify = F)
  res = do.call(rbind, res)
  return(res)
}, simplify = F))
# View(IDS)
write.table(IDS, file = "temp/histoneChipENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)
metafile = do.call(rbind,sapply(tissues, function(tissue){
  temp = read.table(paste0("data/rawdata/histoneENCODE/", tissue, "/metafile.txt"))
  return(c(tissue, temp))
}, simplify = F))
write.table(metafile, file = "temp/histoneChipENCODE_metafile.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = F)
#####

# DNAse ENCODE (=DNAaccessibility)#####
IDS = do.call(rbind, sapply(tissues, function(tissue){
  meta = read.table(paste0("data/rawdata/DNAseENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T, quote = "")
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type %in% c("bigWig", "bed") , ]
  
  filesBW = meta$File.accession[meta$File.type == "bigWig"]
  filesBed = meta$File.accession[meta$File.type == "bed"]
  
  res = rbind(cbind(tissue,  "bigWig", paste(filesBW, collapse = ", ")), 
              cbind(tissue,  "bed", paste(filesBed, collapse = ", ")))
  res = res[res[,3] != "",]
  return(res)
}, simplify = F))
# View(IDS)
write.table(IDS, file = "temp/DNAseENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)
metafile = do.call(rbind,sapply(tissues, function(tissue){
  temp = read.table(paste0("data/rawdata/DNAseENCODE/", tissue, "/metafile.txt"))
  return(c(tissue, temp))
}, simplify = F))
write.table(metafile, file = "temp/DNAseENCODE_metafile.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = F)
#####

# ATACseq ENCODE #####
IDS = do.call(rbind, sapply(tissues, function(tissue){
  if(!file.exists(paste0("data/rawdata/ATACseqENCODE/", tissue, "/metadata.tsv"))){
    return(NULL)
  }
  meta = read.table(paste0("data/rawdata/ATACseqENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[grep(x = meta$Audit.ERROR,pattern = "extreme",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "severe",invert = T),]
  meta = meta[grep(x = meta$Audit.NOT_COMPLIANT,pattern = "extreme",invert = T),]
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$File.type %in% c("bigWig", "bed") , ]
  
  filesBW = meta$File.accession[meta$File.type == "bigWig"]
  filesBed = meta$File.accession[meta$File.type == "bed"]
  
  res = rbind(cbind(tissue,  "bigWig", paste(filesBW, collapse = ", ")), 
              cbind(tissue,  "bed", paste(filesBed, collapse = ", ")))
  res = res[res[,3] != "",]
  return(res)
}, simplify = F))
# View(IDS)
write.table(IDS, file = "temp/ATACseqENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)
metafile = do.call(rbind,sapply(tissues, function(tissue){
  if(!file.exists(paste0("data/rawdata/ATACseqENCODE/", tissue, "/metafile.txt"))){
    return(NULL)
  }
  temp = read.table(paste0("data/rawdata/ATACseqENCODE/", tissue, "/metafile.txt"))
  return(c(tissue, temp))
}, simplify = F))
write.table(metafile, file = "temp/ATACseqENCODE_metafile.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = F)
#####

# HiC #####
IDS = do.call(rbind, sapply(tissues, function(tissue){
  meta = read.table(paste0("data/rawdata/HiCENCODE/", tissue, "/metadata.tsv"),
                    header = T, sep = "\t", as.is = T)
  meta = meta[meta$File.Status == "released",]
  meta = meta[meta$Output.type %in% c("genome compartments",
                                      "mapping quality thresholded contact matrix"),]
  
  filesComp = meta$File.accession[meta$Output.type == "genome compartments"]
  filesMatrix = meta$File.accession[meta$Output.type == "mapping quality thresholded contact matrix"]
  
  res = rbind(cbind(tissue,  "genome compartments", paste(filesComp, collapse = ", ")), 
              cbind(tissue,  "mapping quality thresholded contact matrix", paste(filesMatrix, collapse = ", ")))
  res = res[res[,3] != "",]
  return(res)
}, simplify = F))
# View(IDS)
write.table(IDS, file = "temp/HiCENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)
metafile = do.call(rbind,sapply(tissues, function(tissue){
  temp = read.table(paste0("data/rawdata/HiCENCODE/", tissue, "/metafile.txt"))
  return(c(tissue, temp))
}, simplify = F))
write.table(metafile, file = "temp/HiCENCODE_metafile.tsv", sep = "\t", 
            quote = F, row.names = F, col.names = F)
#####


# Replication ######
IDS = sapply("x", function(dummy){
  files =read.table("data/rawdata/replication/files.txt", sep = "\t")
  temp = lapply(files$V2,function(x){
    x = strsplit(x,split="; ")[[1]]
    x = do.call(rbind,strsplit(x,split="="))
    setNames(x[,2], x[,1])
  })
  tags = unique(unlist(sapply(temp, names)))
  meta = t(sapply(temp, function(x){
    x[tags]
  }))
  colnames(meta) = tags
  meta = data.frame(fileName = files$V1, meta)
  res = meta[meta$view %in% c("WaveSignal", "Peaks", "Valleys"),
             c("fileName", "view", "cell",  
               "dccAccession", "subId",  "type")]
  return(res)
}, simplify = F)[[1]]
write.table(IDS, file = "temp/ReplicationENCODE_IDs.tsv", sep = "\t", quote = F, row.names = F)

#####


# methylation ######
IDS = do.call(rbind,sapply(tissues, function(tissue){
  files = list.files(paste0("data/rawdata/methylationGSE186458/", tissue))
  return(cbind(tissue, paste(files, collapse = ", ")))
}, simplify = F))
write.table(IDS, file = "temp/methylation_IDs.tsv", sep = "\t", quote = F, row.names = F)

######