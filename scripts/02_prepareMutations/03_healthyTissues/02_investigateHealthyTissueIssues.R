# preparation
# library(rtracklayer)
library(GenomicRanges)
chrs = paste0("chr", c(1:22))
chrTranslator = setNames(c(chrs, chrs),
                         c(chrs, 1:22))

tissues2Test = c("prostate", "liver", "ureter", "rectum", "skin", "pancreas", "esophagus")
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # gr
#####

# for the example tissues, load TPs #####
TP_list = sapply(tissues2Test, function(tissue){
  TPs = read.table(paste0("data/rawdata/SomaMutDB/",tissue,".tsv"),
                            sep = "\t", header = T, comment.char = "")
  colnames(TPs)[colnames(TPs) == "X.chr"] = "chr"
  TPs$chr = chrTranslator[TPs$chr]
  TPs = TPs[!is.na(TPs$chr),]
  TPs <- TPs[nchar(TPs$ref) == 1 & nchar(TPs$alt) == 1,]
  posToFilter = paste(TPs$chr, TPs$pos, sep = "_")
  dupl = duplicated(posToFilter) 
  TPs = TPs[!dupl,]
  TPsGR = GRanges(seqnames=TPs$chr, 
                  ranges=IRanges(start=TPs$pos, end=TPs$pos))
  ov = is.na(findOverlaps(query=TPsGR, subject=gr, select = "first"))
  TPs = TPs[ov,c("chr", "pos", "ref", "alt")]
}, simplify = F)
#####

# load UCSC problematic region annotations #####
# downloaded into data/rawdata/genome/UCSCproblematic
# wget ftp://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/problematic/*
# probFiles = list.files("data/rawdata/genome/UCSCproblematic/", pattern = ".bb")
# command = paste0("lib/bigBedToBed -tsv data/rawdata/genome/UCSCproblematic/", probFiles, 
#                  " data/rawdata/genome/UCSCproblematic/",
#                  sapply(strsplit(probFiles, split = ".", fixed = T), function(x){x[1]}), ".bed",
#                  collapse = "; ")
# system(command, intern = T)
probBeds = list.files("data/rawdata/genome/UCSCproblematic/", pattern = ".bed")
problematic = sapply(probBeds, function(fil){
  print(fil)
  probGR = read.table(paste0("data/rawdata/genome/UCSCproblematic/", fil),
                      sep = "\t", header = T, quote = "")
  probGR = probGR[probGR$chrom %in% chrs,]
  return(probGR)
}, simplify = F)
#####


# visualize #####
pdf("fig/healthyTissuesWGS/problematicTissuesMutsAlongGenome.pdf")
par(mfrow = c(length(tissues2Test),1),mar = c(2,2,1,1), oma = c(1,2,3,0))
toPlot = c("comments.bed",  "encBlacklist.bed", 
           "ngsProblemHigh.bed", "sangerDeadZone.bed")
dumpVar = sapply(chrs, function(cr){
  dumpVar2 = sapply(tissues2Test, function(tissue){
    TPs = TP_list[[tissue]]
    temp = hist(TPs$pos[TPs$chr == cr], breaks = 200, main = tissue)
    dumpVar2 = sapply(seq_len(length(toPlot)), function(i){
      track = problematic[[toPlot[i]]]
      track = track[track$chrom == cr,]
      yPos = seq(0,max(temp$counts), length.out = length(toPlot))
      if(nrow(track)>=1){
        segments(x0 = track$chromStart, x1 = track$chromEnd, 
                 y0 = yPos[i], y1 = yPos[i], 
                 lwd = 3, col = i+1)
      }
    })
  })
  mtext(cr, side = 3, outer = T)
  mtext("position", side = 1, outer = T)
  mtext("n Mutations", side = 2, outer = T)
})
dev.off()
# chromosome 2 #
temp = hist(TP_list$ureter[TP_list$ureter$chr == "chr2","pos"], breaks = 300)
temp$breaks[which.max(temp$counts):(which.max(temp$counts)+1)]
# chr2:179000000-180000000
track = problematic$ngsProblemHigh.bed
track[track$chrom == "chr2" & track$chromStart>1.7e8 & track$chromEnd<1.8e8,]
# --> TTN (titin) which lies at chr2:179,390,716-179,672,150
par(mfrow = c(length(tissues2Test),2),mar = c(2,2,1,1), oma = c(1,2,3,0))
dumpVar2 = sapply(tissues2Test, function(tissue){
  TPs = TP_list[[tissue]]
  hist(TPs$pos[TPs$chr == "chr2"], breaks = 200, main = tissue)
  TPs = TPs[TPs$chr == "chr2" & (TPs$pos < 179390716 | TPs$pos > 179672150),]
  hist(TPs$pos, breaks = 200, main = tissue)
})
# yes, it's due to TTN!

# chromosome 5
temp = hist(TP_list$ureter[TP_list$ureter$chr == "chr5","pos"], breaks = 300)
peak = which(temp$counts>50)
temp$breaks[c(peak, peak+1)]
# chr5:140000000-141000000
track = problematic$comments.bed
track[track$chrom == "chr5" & track$chromStart>140500000 & track$chromEnd<141000000,]
# --> PCDHA or PCDHB clusters 
cluster = GRanges(track[track$chrom == "chr5",1:4])
par(mfrow = c(length(tissues2Test),2),mar = c(2,2,1,1), oma = c(1,2,3,0))
dumpVar2 = sapply(tissues2Test, function(tissue){
  TPs = TP_list[[tissue]]
  hist(TPs$pos[TPs$chr == "chr5"], breaks = 200, main = tissue)
  TPsGR = GRanges(seqnames = TPs$chr, ranges = IRanges(TPs$pos, width = 1))
  ov = is.na(findOverlaps(query=TPsGR, subject=cluster, select = "first"))
  TPsfiltered = TPs[ov,]
  hist(TPsfiltered$pos[TPsfiltered$chr == "chr5"], breaks = 200, main = tissue)
})
# yes, it's due to the clusters

# chromosome 7
temp = hist(TP_list$ureter[TP_list$ureter$chr == "chr7","pos"], breaks = 1000,
            xlim = c(8e7, 1.6e8))
dumpVar2 = sapply(seq_len(length(problematic)), function(i){
  track = problematic[[i]]
  track = track[track$chrom == cr,]
  yPos = seq(0,max(temp$counts), length.out = length(problematic))
  if(nrow(track)>=1){
    segments(x0 = track$chromStart, x1 = track$chromEnd, 
             y0 = yPos[i], y1 = yPos[i], 
             lwd = 1, col = i+1)
  }
})
peak = which(temp$counts>20)
temp$breaks[c(peak, peak+1)]
# chr7:100600000-100800000
# chr7:149400000-149600000

#####


# test performance if we remove problematic chromosomes #####
library(ROCR)
perfsAllTissueModelReduced = sapply(tissues2Test, function(tissue){
  # created in scripts/02_prepareMutations/05b_healthyTissues_WGS.R :
  load(paste0("temp/tempRes_", tissue, ".RData")) # res
  res = res[!res$chrom %in% c("chr2", "chr5", "chr6", "chr7", "chr8", "chr9"),]
  temp = prediction(pred = res$pred,labels = res$labels)
  roc = performance(temp,  "tpr", "fpr")
  auc = performance(temp, "auc")
  pr = performance(temp, "prec", "rec")
  return(list(roc = roc, pr = pr, auc = auc))
}, simplify = F)
load("data/Modeling/healthyTissues/perfsAllTissueModel_WGS.RData") # perfsAllTissueModel

png("fig/healthyTissuesWGS/problematicTissuesPerfWithoutProbChroms.png", 
    height = 1200, width = 400, pointsize = 12)
par(mfrow = c(7,2), mar = c(1,1,2,1), oma = c(3,3,0,0))
dumpVar = sapply(tissues2Test, function(tissue){
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "FPR", ylab = "TPR", las = 1, main = tissue)
  abline(0,1, col = "dark grey", lty = 2, lwd = 1.5)
  # general model
  lines(perfsAllTissueModel[[tissue]]$roc@x.values[[1]],
        perfsAllTissueModel[[tissue]]$roc@y.values[[1]])
  lines(perfsAllTissueModelReduced[[tissue]]$roc@x.values[[1]],
        perfsAllTissueModelReduced[[tissue]]$roc@y.values[[1]], 
        col = "red", lty = 2)
  
  plot(NA, xlim = c(0,1), ylim = c(0,1),
       xlab = "Recall", ylab = "Precision", las = 1,  main = tissue)
  lines(perfsAllTissueModel[[tissue]]$pr@x.values[[1]],
        perfsAllTissueModel[[tissue]]$pr@y.values[[1]])
  lines(perfsAllTissueModelReduced[[tissue]]$pr@x.values[[1]],
        perfsAllTissueModelReduced[[tissue]]$pr@y.values[[1]], 
        col = "red", lty = 2)
})
dev.off()



# Check again after re-processing of mutations: #####
tissues2Test = c("Prostate", "Liver", "Ureter",  "Skin", "Pancreas", "Esophagus") # "Rectum" is no longer in WGS set
chrs = paste0("chr", 1:22)
load("data/processedData/chrLengths.RData")
TP_list = sapply(tissues2Test, function(tissue){
  load(paste0("data/MutTables/SomamutDB/", tissue, "_WGS.RData")) # Muts
  return(Muts)
}, simplify = F)
pdf("fig/healthyTissuesWGS/checkBiasesMutsAlongGenome.pdf")
par(mfrow = c(length(tissues2Test),1),mar = c(2,2,1,1), oma = c(1,2,3,0))
dumpVar = sapply(chrs, function(cr){
  dumpVar2 = sapply(tissues2Test, function(tissue){
    TPs = TP_list[[tissue]]
    temp = hist(TPs$pos[TPs$chr == cr & TPs$mutated == "1"], breaks = 200, 
                main = tissue, xlim = c(0,chrLengths[cr]))
  })
  mtext(cr, side = 3, outer = T)
  mtext("position", side = 1, outer = T)
  mtext("n Mutations", side = 2, outer = T)
  return(NULL)
})
dev.off()
#####