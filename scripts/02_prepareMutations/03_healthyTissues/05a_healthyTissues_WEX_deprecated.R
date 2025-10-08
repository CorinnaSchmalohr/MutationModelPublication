### Packages ######
library(GenomicRanges)
library(stringr)

# library(parallel)
###############

# preparation #####
# nThread = 10
dir.create("data/MutTables/healthyTissuesWEX/", showWarnings = F)
load("data/predictors/wgEncodeDacMapabilityConsensusExcludable.RData") # gr
load("data/processedData/codExons_noUTR_filtered.RData") # codExons_filtered

tissues =c("adipocytes", "adrenal_gland", "bladder", 
           "bonemarrow", "brain", "breast", "colon", 
           "esophagus", "fibroblast", "heart",  "kidney",
           "liver", "lung", "pancreas", "placenta", "prostate", "rectum", 
           "skeletal_muscle", "skin", "small_intestine", "spleen", 
           "stomach", "testicle", "thyroid", "tonsil", "ureter") #"iPSC","blood", "embryonic_stem_cell", "endometrium"(has only indels)
chrs = paste0("chr", c(1:22))
chrTranslator = setNames(c(chrs, chrs),
                         c(chrs, 1:22))
load("data/processedData/pentamers.RData")
rownames(fivemers) = pent2context[fivemers[,1]]
# load sequences for TN generation
exonSeqs = read.table("data/processedData/codExons_noUTR_filtered.withsequences.bed")
# in this file, covered regions are extended by two bases in 
# both direction, in order to be able to find the 5mers.
colnames(exonSeqs) = c("chr", "start", "end", "id","sequence")
exonSeqs = exonSeqs[exonSeqs$chr %in% chrs,]


#####

# iterate through tissues and prepare data #####
sapply(tissues, function(tissue){
  print(tissue)
  TPs = read.table(paste0("data/rawdata/SomaMutDB/",tissue,".tsv"),
                   sep = "\t", header = T, comment.char = "")
  
  # Remove indels 
  TPs <- TPs[nchar(TPs$ref) == 1 & nchar(TPs$alt) == 1,]
  
  # only chrs 1-22
  colnames(TPs)[colnames(TPs) == "X.chr"] = "chr"
  TPs$chr = chrTranslator[TPs$chr]
  TPs = TPs[!is.na(TPs$chr),]
  
  # subset to exome regions
  codExons_filtered = GRanges(codExons_filtered)
  TPsGR = GRanges(seqnames= TPs$chr, 
                  ranges=IRanges(start = TPs$pos,
                                 end=TPs$pos))
  ov = !is.na(findOverlaps(query=TPsGR, subject=codExons_filtered, select = "first"))
  # print(table(ov))
  TPs = TPs[ov,]
  
  
  # store mutations to exclude during TN generation
  posToFilter = paste(TPs$chr, TPs$pos, sep = "_")
  
  # exclude positions that were mutated more than once
  dupl = duplicated(posToFilter) 
  TPs = TPs[!dupl,]
  
  # get sequence context for positions
  TPs_bed = data.frame(TPs$chr, TPs$pos-3, TPs$pos+2)
  TPs_bed = format(TPs_bed, trim = T)
  tmpfile = tempfile(pattern = tissue, fileext = ".bed", tmpdir = "temp")
  write.table(TPs_bed, file=tmpfile, quote=F, col.names=F, row.names=F, sep="\t")
  cmd = paste0("bedtools getfasta -fi data/rawdata/GRCh37.primary_assembly.genome.fa ", 
               "-bed  ",tmpfile,
               " -name -tab | cut -f 2")
  context = system(cmd, intern = T)
  TPs$context = context
  if(!all.equal(substring(context, 3,3), TPs$ref)){
    error("context is not matching reference base. check script")
  }
  file.remove(tmpfile)
  
  #  filter for callable regions
  # this step is actually unnecessary since the codExons_filtered is already filtered the same way. Leaving in for consistency
  TPsGR = GRanges(seqnames=TPs$chr, 
                  ranges=IRanges(start=TPs$pos, end=TPs$pos))
  ov = is.na(findOverlaps(query=TPsGR, subject=gr, select = "first"))
  TPs = TPs[ov,c("chr", "pos", "ref", "alt","context")]
  if(nrow(TPs)==0){
    print("no mutations left for this tissue")
    return(NULL)
  }
  # Get TNs
  TNs = lapply(unique(TPs$chr), function(chr){ 
    # subset data to chromosome
    TPs_chr = TPs[TPs$chr == chr,]
    exonSeqs_chr = exonSeqs[exonSeqs$chr == chr,]
    pents = TPs_chr$context
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
  TNs = do.call(rbind,TNs)
  TNs = as.data.frame(TNs)
  TNs$pos = as.integer(TNs$pos)
  colnames(TNs) = c("chr", "pos", "ref", "context")
  
  # save
  TNs$alt = NA
  TNs$mutated = 0
  TPs$mutated = 1
  Muts = rbind(TPs, TNs)
  Muts = Muts[order(Muts[,1], Muts[,2]),]
  Muts = Muts[,c("chr", "pos", "ref", "alt", "context", "mutated")]
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                  ids)
  save(Muts, file = paste0("data/MutTables/healthyTissuesWEX/", tissue, "_SomaMutDB.RData"))
  
  ids = do.call(paste,c(Muts, sep = "_"))
  MutsBed = cbind(Muts[,1], as.integer(Muts[,2]-1), as.integer(Muts[,2]), 
                  ids)
  write.table(MutsBed, file = paste0("data/MutTables/healthyTissuesWEX/", tissue, "_SomaMutDB.bed"),
              col.names = F, row.names = F, sep = "\t", quote = F)
})
#####
