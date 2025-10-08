# preparation #####
tissue2Cancer = list("lung" = "Lung adenocarcinoma",
                     "breast" = "Breast invasive carcinoma",
                     "skin" = "Skin Cutaneous Melanoma",
                     "colon" = c("Colon adenocarcinoma", "Rectum adenocarcinoma"),
                     "ovary" = "Ovarian serous cystadenocarcinoma",
                     "kidney" = "Kidney renal clear cell carcinoma",
                     "prostate" = "Prostate adenocarcinoma",
                     "esophagus" = "Esophageal carcinoma",
                     "liver" = "Liver hepatocellular carcinoma",
                     "brain" = "Brain Lower Grade Glioma")
chrs = paste0("chr", c(1:22))
library(stringr)
library(GenomicRanges)
library(Biostrings)
#####

# load the pancan mutations #####
muts = read.table("data/rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz",
                  header = T, na.strings=".", quote="", as.is=T)
# get  cancer IDs
cancerIDs = read.table("data/rawdata/pancan/TCGA-CDR-SupplementalTableS1_TSScode_info.csv",
                       sep = ",", header = T, as.is = T, skip=2, na.strings="")
cancerID2cancer = setNames(cancerIDs$Disease.type, cancerIDs$TSS.Code)
cancerCode = substr(muts$Tumor_Sample_Barcode, 6,7)
muts$cancerType = cancerID2cancer[cancerCode]
muts = muts[muts$cancerType %in% unlist(tissue2Cancer) & muts$Variant_Type == "SNP",]
muts$Chromosome = paste0("chr", muts$Chromosome)
muts = muts[muts$Chromosome %in% chrs,]
#####

# extract positions to exclude during TN sampling #####
posToFilterTotal = sapply(names(tissue2Cancer), function(tissue){
  temp = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
  cbind(temp$Chromosome,
        paste(temp$Chromosome, temp$Start_Position, sep = "_"))
})
#####


# filter out variants with too little depth #####
muts = muts[muts$n_depth>=8 & muts$t_depth >= 14,]
# these cutoffs are based on what was used for the "sufficient coverage" files.
#####


# filter out variants with >1% population frequency  #####
mafs = muts[,c("GMAF", "AFR_MAF", "AMR_MAF", "EAS_MAF",
               "EUR_MAF", "SAS_MAF", "AA_MAF", "EA_MAF")]
mafsFilter = apply(mafs,1, function(x){
   y = gsub("[T,C,G,A,-]:", replacement="",x=x)
   z = as.numeric(unlist(strsplit(split=",", x=y)))
   if(all(is.na(z))){
      return(T)
   }  else if(max(z, na.rm = T) < 0.01){
      return(T)
   } else{
      return(F)
   }
})
muts = muts[mafsFilter,]
#####


# subset to exome regions #####
load("data/processedData/codExons_noUTR_filtered.RData")
codExons_filtered = GRanges(codExons_filtered)
mutGR = GRanges(seqnames= muts$Chromosome, 
                ranges=IRanges(start = muts$Start_Position,
                               end=muts$End_Position))
ov = !is.na(findOverlaps(query=mutGR, subject=codExons_filtered, select = "first"))
print(table(ov))
muts = muts[ov,]
# save to use for downsampling
save(muts, file = "data/MutTables/exomeTrainData/muts.RData")
#####


# load sequences for TN generation #####
exonSeqs = read.table("data/processedData/codExons_noUTR_filtered.withsequences.bed")

# in this file, covered regions are extended by two bases in 
# both direction, in order to be able to find the 5mers.
colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]
load("data/processedData/pentamers.RData")
rownames(fivemers) = fivemers[,1]
#####



# for each tissue #####
set.seed(235)
for(tissue in names(tissue2Cancer)){
   print(tissue)
   # subset to tissue
   submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
   
   # exclude positions that were mutated more than once, 
   # to exclude possible selection
   positions = paste(submuts$Chromosome, submuts$Start_Position, sep = "_")
   multiple = duplicated(positions) | duplicated(positions, fromLast=T)
   submuts = submuts[!multiple,]
   
   # get TNs #####
   TNs = lapply(unique(submuts$Chromosome), function(chr){ 
     cat(chr, ' ')
     # subset data to chromosome
     TPs_chr = submuts[submuts$Chromosome == chr,]
     exonSeqs_chr = exonSeqs[exonSeqs$chr == chr,]
     context = TPs_chr$CONTEXT
     contextL = unique(str_length(context))
     pents = substr(context,start=(contextL-1)/2-1, stop=(contextL-1)/2+3)
     posToFilter = posToFilterTotal[[tissue]][posToFilterTotal[[tissue]][,1] == chr,]
     TNs_chr = lapply(unique(pent2context[pents]), function(pent){
       # cat(pent, ' ')
       pent2search = fivemers[pent,]
       nmatch1 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[1])
       nmatch2 = str_count(string = exonSeqs_chr$sequence, pattern = pent2search[2])
       nmatch = nmatch1 + nmatch2
       weights = nmatch/sum(nmatch)
       toRm = rep(1,sum(pents %in% pent2search)) 
       pos = NULL
       while(length(toRm) > 0){
         TNs_sample = sample(1:nrow(exonSeqs_chr),
                             size = length(toRm),
                             prob = weights, replace = T)
         new = t(sapply(TNs_sample, function(i){
           # cat(i,' ')
           gene = unlist(exonSeqs_chr[i,])
           matches1 = cbind(gregexpr(pattern = pent2search[1], text = gene["sequence"])[[1]],pent2search[1])
           matches2 = cbind(gregexpr(pattern = pent2search[2], text = gene["sequence"])[[1]],pent2search[2])
           matches = rbind(matches1, matches2)
           matches = matches[matches[,1] != "-1",, drop = F]
           if(nrow(matches) == 1){
             sampMatch = matches[1,]
           } else{
             sampMatch = matches[sample(1:nrow(matches),size = 1),]
           }
           samp = as.integer(sampMatch[1])
           c(gene["chr"], 
             pos = (as.numeric(gene["start"]) + samp +2),
             ref = substr(gene["sequence"],start = samp+2, stop = samp+2),
             context = sampMatch[2])
         }))
         pos = rbind(pos,new)
         temp = paste(pos[,"chr"], pos[,"pos"], sep = "_")
         #  make sure that none of the positions are actually overlapping with existing positions:
         toRm = which(temp %in% posToFilter | duplicated(temp))
         if(length(toRm) > 0){
           pos = pos[-toRm,]
         }
       }
       return(pos)
     })
     TNs_chr = do.call(rbind,TNs_chr)
     
   })
   TNs = do.call(rbind, TNs)
   TNs = data.frame(chr=TNs[,"chr"], pos = as.integer(TNs[,"pos"]),
                                     ref=TNs[,"ref.sequence"], alt = NA, context=TNs[,"context"],
                                     stringsAsFactors = F)

   # combine TPs and TNs
   context = submuts$CONTEXT
   contextL = unique(str_length(context))
   exomemuts = cbind(submuts[,c("Chromosome", "Start_Position", 
                                "Reference_Allele","Tumor_Seq_Allele2")], 
                     substr(context,start=(contextL-1)/2-1, stop=(contextL-1)/2+3))
   colnames(exomemuts) = colnames(TNs)
   TNs$mutated = 0
   exomemuts$mutated = 1
   Muts = rbind(TNs,exomemuts)
   Muts = Muts[order(Muts[,1], Muts[,2]),]
   
   
   # save
   ids = do.call(paste,c(Muts, sep = "_"))
   MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                   ids)
   save(Muts, file = paste0("data/MutTables/exomeTrainData/", tissue, "_Muts.RData"))
   write.table(MutsBed, file = paste0("data/MutTables/exomeTrainData/", tissue, "_Muts.bed"),
               col.names = F, row.names = F, sep = "\t", quote = F)
}
#####

# get number of mutations per sample per tissue #####
nPerTissue = sapply(names(tissue2Cancer), function(tissue){
   print(tissue)
   # subset to tissue
   submuts = muts[muts$cancerType %in% tissue2Cancer[[tissue]],]
   table(submuts$Tumor_Sample_Barcode)
})
nPerTissue = nPerTissue[sort(names(nPerTissue))]
save(nPerTissue, file="data/MutTables/exomeTrainData/nPerTissue.RData")
####
