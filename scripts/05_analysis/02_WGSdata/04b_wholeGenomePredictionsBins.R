
library(GenomicRanges)
methods = c("RFpreds", "contextPreds", "mult", "combOdds")
source("scripts/05_analysis/00_NamesAndColors.R")

# define bins #####
load("/cellfile/datapublic/ypaul1/Mutations/results/GenomePredictions/autosomes_non_excludable_regions.RData") # result_df
regionsGR =GRanges(result_df)
regionsGR <- reduce(regionsGR)
windowSizes = c("1Mb" = 1000000, "100kb" = 100000, "1kb"  =1000, "100bp" = 100)
windows = sapply(names(windowSizes), function(size){
  bins = slidingWindows(regionsGR, width = windowSizes[size], step =  windowSizes[size]) 
  bins = unlist(bins)
  bins = bins[width(bins) == windowSizes[size]]
  bins = bins[seqnames(bins) == "chr22"]
  return(bins)
}, simplify = F)
#####
windowSizes = c(windowSizes, "1bp" = 1)


# iterate through windowsizes
# for each, average predictions, sum up mutations counts, then compute performance
perf = sapply(names(windowSizes), function(size){
  print(size)
  perfTissues = sapply(WGStissues, function(tissue){
    print(tissue)
    if(tissue == "ovary"){
      mutationCols = "cancerMuts"
    } else{
      mutationCols = c("cancerMuts","healthyMuts")
    }
    # load exome-wide prediction table
    load(paste0("data/Modeling/WholeGenomeData/combinedPredictions/combinedPredictions_", 
                tissue, "_chr","22",".RData")) #data
    if(size == "1bp"){
      tempDat = na.omit(mcols(data))
      scores = as.matrix(tempDat[,methods])
      mutations = as.matrix(tempDat[,mutationCols])
      res = list(correlation = cor(scores, mutations, method = "spearman"),
                 correlation_pval = apply(mutations,2,function(x){
                   apply(scores,2,function(y){
                     cor.test(x,y, method = "spearman", exact = F)$p.value
                   })
                 }))
      return(res)
    }
      
    # assign positions to bins
    hits  =  findOverlaps(query = windows[[size]], subject = data)
    bin_idx <- queryHits(hits)
    # compute mean score for each bin
    avgScores = sapply(methods, function(varname){
      scores <- mcols(data)[,varname][subjectHits(hits)]
      # Compute mean score per bin using tapply
      meanPerBin <- tapply(scores, bin_idx, mean, na.rm = T)
      return(meanPerBin)
    })
    
    # compute number of mutations for each bin

    binnedMutations = sapply(mutationCols, function(varname){
      scores <- mcols(data)[,varname][subjectHits(hits)]
      # Compute sum of mutations per bin using tapply
      meanPerBin <- tapply(scores, bin_idx, sum)
      return(meanPerBin)
    })
    res = list(correlation = cor(avgScores, binnedMutations, method = "spearman"),
                correlation_pval = apply(binnedMutations,2,function(x){
                  apply(avgScores,2,function(y){
                    cor.test(x,y, method = "spearman", exact = F)$p.value
                  })
                }))
    return(res)
  }, simplify = F)
  return(perfTissues)
}, simplify = F)
save(perf, file = "data/Modeling/WholeGenomeData/perf_Allmethods_binned.RData")



# visualize #####
load(paste0("data/Modeling/WholeGenomeData/perf_Allmethods_", "allTissuesTogether", ".RData")) # perf
singleBP = perf
load("data/Modeling/WholeGenomeData/perf_Allmethods_binned.RData") # perf
plotPerf = function(x, tissue, names.arg = NA,p){
  dens = ifelse(p<=0.05, yes = -1, no = 20)
  temp = barplot(x, density = dens,#ylim = c(0,max(x)),
          las = 1, names.arg = names.arg, col = tissueCols[tissue], xpd = NA) #, horiz = T
  abline(h=0)
  return(temp)
}
png("fig/wholeGenomePredictions/perf_cancer_binned_overview.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(8,5), mar = c(2,2.5,1,1), oma = c(4,8,2,0))
dumpVar = sapply(WGStissues, function(tissue){
  print(tissue)
  dumpVar2 = sapply(rev(names(windowSizes)), function(size){
    if(tissue == "ovary"){
      temp = plotPerf(perf[[size]][[tissue]]$correlation[,1], 
                      tissue = tissue,
                      p = perf[[size]][[tissue]]$correlation_pval[,1])
    }else{
      temp = plotPerf(perf[[size]][[tissue]]$correlation[,"cancerMuts"], 
                      tissue = tissue,
                      p = perf[[size]][[tissue]]$correlation_pval[,"cancerMuts"])
    }
    if(tissue == WGStissues[1])
      title(main=size, line = 1, xpd = NA)
    if(tissue == tail(WGStissues,1))
      text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), 
           labels  = c("RF", "Context", "Mult", "Comb. Odds"), 
           srt = 45, xpd = NA, adj = 1, offset = 5)
    if(size == names(rev(windowSizes))[1]){
      if(tissue == "ovary"){
        text(x = -3.5, y = max(perf[[size]][[tissue]]$correlation[,1])/2,
             labels = t2T[tissue], adj = 1, xpd = NA)
      }else{
        text(x = -3.5, y = max(perf[[size]][[tissue]]$correlation[,"cancerMuts"])/2,
             labels = t2T[tissue], adj = 1, xpd = NA)
      }
    }
      
  })
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 6.6)
dev.off()



png("fig/wholeGenomePredictions/perf_healthy_binned_overview.png", 
    width = 1200, height = 1200, pointsize = 30)
par(mfrow = c(7,4), mar = c(2,2.5,1,1), oma = c(4,8,2,0))
dumpVar = sapply(WGStissues[WGStissues != "ovary"], function(tissue){
  print(tissue)
  dumpVar2 = sapply(rev(names(windowSizes)), function(size){
    temp = plotPerf(perf[[size]][[tissue]]$correlation[,"healthyMuts"], 
                    tissue = tissue,
                    p = perf[[size]][[tissue]]$correlation_pval[,"healthyMuts"])
    if(tissue == WGStissues[1])
      title(main=size, line = 1, xpd = NA)
    if(tissue == tail(WGStissues,1))
      text(x = temp[,1], y = par("usr")[3]-(par("usr")[4]*0.2), 
           labels  = c("RF", "Context", "Mult", "Comb. Odds"), 
           srt = 45, xpd = NA, adj = 1, offset = 5)
    if(size == names(rev(windowSizes))[1])
      text(x = -2.8, y = median(axTicks(2)),
           labels = t2T[tissue], adj = 1, xpd = NA)
  })
})
mtext("Correlation with binned mutation count", side = 2, outer = T, line = 6.6)
dev.off()
######

# compute R squared ######
Rsq = sapply(names(windowSizes), function(size){
  print(size)
  perfTissues = sapply(WGStissues, function(tissue){
    print(tissue)
    if(tissue == "ovary"){
      mutationCols = "cancerMuts"
    } else{
      mutationCols = c("cancerMuts","healthyMuts")
    }
    # load exome-wide prediction table
    load(paste0("data/Modeling/WholeGenomeData/combinedPredictions/combinedPredictions_", 
                tissue, "_chr","22",".RData")) #data
    if(size == "1bp"){
      tempDat = na.omit(mcols(data))
      scores = as.matrix(tempDat[,methods])
      mutations = as.matrix(tempDat[,mutationCols])
      res = apply(scores,2,function(y){
        cor(mutations[,"cancerMuts"], y)^2 
      })
      return(res)
    }

    # assign positions to bins
    hits  =  findOverlaps(query = windows[[size]], subject = data)
    bin_idx <- queryHits(hits)
    # compute mean score for each bin
    avgScores = sapply(methods, function(varname){
      scores <- mcols(data)[,varname][subjectHits(hits)]
      # Compute mean score per bin using tapply
      meanPerBin <- tapply(scores, bin_idx, mean, na.rm = T)
      return(meanPerBin)
    })
    
    # compute number of mutations for each bin
    binnedMutations = sapply(mutationCols, function(varname){
      scores <- mcols(data)[,varname][subjectHits(hits)]
      # Compute sum of mutations per bin using tapply
      meanPerBin <- tapply(scores, bin_idx, sum)
      return(meanPerBin)
    })
    
    # compute Rsquared for cancer mutations only
    res = apply(avgScores,2,function(y){
      cor(binnedMutations[,"cancerMuts"], y)^2 
    })
    return(res)
  }, simplify = F)
  return(perfTissues)
}, simplify = F)
save(Rsq, file = "data/Modeling/WholeGenomeData/perfRsq_Allmethods_binned.RData")

######