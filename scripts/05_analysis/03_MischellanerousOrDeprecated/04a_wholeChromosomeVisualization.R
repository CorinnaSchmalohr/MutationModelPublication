
# get map of i to positions #####
# iPos = sapply(1:32, function(i){
#    cat(i,' ')
#    load(paste0("data/rdata/exomeMuts/exomeMuts_part",i,".RData"))
#    temp = split(subMuts, subMuts$chr)
#    do.call(rbind,lapply(names(temp), function(cr){
#       data.frame(chr = cr,
#                  start = min(temp[[cr]]$pos),
#                  end = max(temp[[cr]]$pos))
#    }))
# }, simplify = F); cat('\n')
# names(iPos) = 1:32
# save(iPos, file = "data/rdata/iPos.RData")
######



# preparation #####
library(readxl)
dir.create("fig/exomeVisualization", showWarnings=F)
load("data/rdata/iPos.RData")
gtf = read.table("data/rawdata/gencode.v43lift37.basic.annotation.gtf.gz", sep = "\t")
regions = read.table("data/rawdata/WEXregions.csv", sep = ",", header=T)
tissues = c("brain","breast", "colon","esophagus","kidney", "liver", "lung",
            "ovary", "prostate", "skin")
tissueCols = setNames(rainbow(length(tissues)), tissues)
baseCols = c("A" = "orange", "C"="blue", "G" = "green", "T" = "red")
darkCol = function(x,n){ # function to create a darkened version of a color
   sapply(x, function(y){
      cl = colorRampPalette(c(y, "black"))(n+1)[1:n]
   })
}
lightCol = function(x,alpha){
   x = col2rgb(x)
   rgb(red=x[1,], green=x[2,], blue=x[3,], alpha=alpha, maxColorValue=255)
}
#####


# get list of variables, grouped by variable type #####
predictors = unique(unlist(sapply(tissues, function(tissue){
   tab = read_xlsx("data/rawdata/dataMapping_old.xlsx", 
                   sheet=tissue, col_names=T)
   tab = as.data.frame(tab)
   return(tab$abbreviation)
})))
predictors = predictors[!predictors %in% c("ConsensusExcludable", 
                                           "repeatMasker",
                                           "tandemRepeatFinder", "context")]
predictors = predictors[!predictors %in% c("ZBTB33_100bp", "YY1_100bp", "TAF1_100bp", "SP1_100bp",
                                           "RXRA_100bp", "REST_100bp", "RAD21_100bp",
                                           "NR2F2_100bp", "MAX_100bp", "JUND_100bp", "HNF4G_100bp",
                                           "HNF4A_100bp", "GABPA_100bp", "FOXA2_100bp", 
                                           "FOXA1_100bp", "EGR1_100bp", "ATF3_100bp")]
binVars = c("Replication_valleys", "Replication_peaks", 
            "direct_repeats", "g_quadruplex", "aPhased_repeats",
            "inverted_repeats", "mirror_repeats", 
            "short_tandem_repeats", "zDNA")
catVars = c("context_precedingBase", "context_ref", "context_followingBase")
numVars = predictors[!predictors %in% c(binVars, catVars)]
numVarsGrouped = sapply(c("GCcontent", "dist", "Replication", "HiC",
                          "DNAse", "methylation","H3K27me3","H3K27ac",
                          "H3K36me3", "H3K4me1", "H3K4me3", "H3K9me3", 
                          "H3K9ac", "expression","TF", "eQTL", "conservation", 
                          "map", "effect", "CTCF",  "PolR2A", "EP300"), function(x){
   numVars[grep(x,numVars, ignore.case=T)] 
})
names(numVarsGrouped)[names(numVarsGrouped) == "map"] = "mappability"
predicCols = setNames(rainbow(length(numVarsGrouped)), names(numVarsGrouped))
numVarsTissue = numVarsGrouped[!names(numVarsGrouped) %in% 
                                  c("GCcontent", "dist", "Replication", 
                                    "TF", "conservation", "mappability", "effect")]
numVarsTissue[["DNAse"]] = c("DNAse_10kb",  "DNAse_100bp")
if( any(!numVars %in% unlist(numVarsGrouped))){error("not all predictors represented")}
predVars = c("RFpreds", "RFcorr","contextPreds", 
             "contextCorr", "combOdds", "combOddsCorr")
predCols = setNames(rainbow(length(predVars)),predVars)
######


# iterate through regions #####
dumpVar = apply(regions,1,function(r){
   print(r[["Gene"]])
   st = as.integer(r[["plotStart"]])
   end = as.integer(r[["plotEnd"]])
   by = as.integer(r[["by"]])
   
   # find out which index to use for data loading
   is = which(sapply(iPos, function(x){
      any(x$chr == r[["chr"]] & 
         ((st >= x$start & st <= x$end) | 
         (end >= x$start & end <= x$end) |
         (st <= x$start & end >= x$end)))
   }))
   
  
   
   # iterate through tissues and load data
   plotData = sapply(tissues, function(tissue){
      cat(tissue, ' ')
      # load predictions
      load(paste0("data/procData/exomeMuts/exomePredictions_",
                                 tissue, ".RData")) # preds
      # scale predictions from 0 to 1 (for easier plotting)
      temp = sapply(colnames(preds[,-(1:3)]),function(x){
         range(preds[,x])
      }) 
      preds = preds[preds$chr == r[["chr"]] & 
                       preds$pos >= st &  
                       preds$pos <= end,]
      preds[,-(1:3)] = sapply(colnames(temp), function(x){
         (preds[,x]-temp[1,x])/(temp[2,x]-temp[1,x])
      })
      # load predictors
      predics = lapply(is, function(i){
         load(paste0("data/procData/exomeMuts/exomeMuts_part",i,"_",
                     tissue,"_processed.RData")) # dat and datpos
         cbind(datpos[,1:2], dat)[datpos$chr == r[["chr"]] & 
                                     datpos$pos >= st &  
                                     datpos$pos <= end,]
      })
      predics = do.call(rbind, predics)
      # numeric predictors are all scaled from 0 to log(2) - rescale from 0 to 1
      predics[,numVars[numVars %in% colnames(predics)]] = 
         sapply(numVars[numVars %in% colnames(predics)], function(x){
         predics[,x] / log(2)
      })
      
      
      # combine predictions and predictors
      pos1 = paste(preds$chr, preds$pos, sep = "_")
      pos2 = paste(predics$chr, predics$pos, sep = "_")
      if(!all.equal(pos1, pos2)){error("non-matching positions")}
      predics$chr = NULL
      predics$pos = NULL
      predics$mutated = NULL
      data = cbind(preds, predics)
      data = data[order(data$pos),]
      
      return(data)
   }, simplify=F); cat('\n')
   
   
   # create plot index for each position so that exons are directly next to
   # each other, with a "buffer" of 50
   pos = plotData[[1]]$pos
   bounds = which(pos[2:length(pos)] - pos[1:(length(pos)-1)] >= 2)
   holes = cbind(holestarts = pos[bounds]+1, holeend = pos[bounds+1]-1)
   exons = cbind(exonstarts = c(min(pos), pos[bounds+1]),
                 exonends = c(pos[bounds], pos[length(pos)]))
   exonL = exons[,2]-exons[,1]+1
   temp = sapply(1:nrow(exons), function(i){
      if(i == 1){
         plotpos = 1:exonL[i]
      } else{
         plotpos = 1:exonL[i] + sum(exonL[seq_len(i-1)]) + 50*(i-1) # buffer of 50
      }
      realpos = exons[i,1]:exons[i,2]
      return(cbind(plotpos, realpos))
   })
   temp = do.call(rbind,temp)
   # translation of real position to index and vice versa
   r2i = approxfun(temp[,2], temp[,1])
   i2r = approxfun(temp[,1], temp[,2])
   rm(bounds,  exonL, temp); gc()
   
   # get annotation for this gene
   if(!is.na(r[["transcriptID"]])){
      annot = gtf[grep(r[["transcriptID"]], gtf$V9),]
      annot = sapply(unique(annot$V3), function(x){
         annot[annot$V3 == x, 4:5]
      }, simplify=F)
   } 
   
   
   # iterate through tissues and create overview figure
   print("overview figure")
   dumpVar2 = sapply(tissues, function(tissue){
      cat(tissue, ' ')
      png(paste0("fig/exomeVisualization/", r[["Gene"]], "_", tissue, ".png"),
          height=2000, width=1400, res = 200)
      {
         # use custom plot margins to adapt dimensions of individual plots.
         par(mar = c(0,9,0,5.7), oma = c(4,0.5,3,0.5))
         layout(rbind(1,2,2,2,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4))
         # gene annotation
         plot(NA, xlim = range(r2i(pos)), ylim = c(-3,2),bty = "n",
              yaxt = "n",xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
         if(!is.na(r[["transcriptID"]])){
            abline(h=1)
            segments(x0=r2i(annot[["transcript"]]$V4), 
                     x1 = r2i(annot[["transcript"]]$V5), 
                     y0=1, lwd = 6, col = "black", lend=1)
            segments(x0=r2i(annot[["exon"]]$V4),
                     x1 = r2i(annot[["exon"]]$V5), 
                     y0=1, lwd = 10, col = "black", lend=1)
            segments(x0=r2i(annot[["UTR"]]$V4), x1 = r2i(annot[["UTR"]]$V5), 
                     y0=1, lwd = 10, col = "blue", lend=1)
            segments(x0=r2i(annot[["CDS"]]$V4), x1 = r2i(annot[["CDS"]]$V5), 
                     y0=1, lwd = 10, col = "red", lend=1)
            if("stop_codon" %in% names(annot)){
               points(x=r2i(annot[["stop_codon"]]$V4), 
                      y=rep(1, nrow(annot[["stop_codon"]])), 
                      pch = "|", col = "red", cex = 2, xpd = NA)
            }
            if("start_codon" %in% names(annot)){
               points(x=r2i(annot[["start_codon"]]$V4), xpd = NA,
                      y=rep(1, nrow(annot[["start_codon"]])), 
                      pch = "|", col = "blue", cex = 1.5)
            }
            legend(x = r2i(max(pos)), y = 1.5, xpd = NA,ncol = 4,
                   legend = c("UTR", "exon", "start", "stop"), xjust = 1,yjust = 0,
                   y.intersp=0.9,bty = "n",
                   pch = c(NA, NA, "|", "|"), col = c(NA, NA, "blue", "red"),
                   fill = c("blue", "red", NA, NA), border=NA)
            legend(x = r2i(min(pos)), y = 1.5, xpd = NA,
                   xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",
                   legend = paste(ifelse(r[["strand"]]=="+", yes="---->", no = "<----"),
                                  "reading direction"),border=NA)
         }
         
         muts = plotData[[tissue]]$mutated
         abline(h=-1.5, col = "grey")
         if(any(muts)){
            segments(x0 = r2i(pos[muts]),y0=-2, y1 = -1)
         }
         text(x = min(r2i(pos)), y = -1.5, 
              labels = "mutations ", 
              xpd = NA, adj = c(1, 0.5))
         
         # predictions
         plot(NA,xlim = range(r2i(pos)),
              ylim = c(1,length(predVars)+1), xaxt ="n", yaxt = "n",
              xlab = "", ylab = "", xaxs = "i", yaxs = "i")
         text(x = min(r2i(pos)), y = 1:length(predVars)+0.5, 
              labels = paste0(predVars," "), 
              xpd = NA, adj = c(1, 0.5), col = darkCol(predCols[predVars], 2)[2,])
         dumpVar3 = sapply(1:length(predVars), function(i){
            temp = plotData[[tissue]][,predVars[i]]
            temp = temp*0.95+0.025
            x = r2i(c(pos[1], pos, pos[length(pos)]))
            y = c(i,temp+i,i)
            polygon(x,y, border=NA, col = predCols[predVars[i]])
            lines(r2i(pos), lwd = 1,
                  temp+i, col = predCols[predVars[i]])
         })
         rect(ybottom=par("usr")[3], ytop=par("usr")[4],
              xleft=r2i(holes[,1]), border=NA,
              xright=r2i(holes[,2]), col="grey")
         rect(ybottom=par("usr")[3], ytop=par("usr")[4],
              xleft=r2i(min(pos)), border=NA,
              xright=r2i(r[["start"]]), col=lightCol("grey", alpha = 90))
         rect(ybottom=par("usr")[3], ytop=par("usr")[4],
              xleft=r2i(r[["end"]]), border=NA,
              xright=r2i(max(pos)), col=lightCol("grey", alpha = 90))
         abline(h=1:length(predVars)+1)
         box()
         
         # 0/1 predictors and sequence context
         tempVars = binVars[binVars %in% colnames(plotData[[tissue]])]
         plot(NA,xlim = range(r2i(pos)),
              ylim = c(0.5,length(tempVars)+length(catVars)+0.5),
              xaxt ="n", yaxt = "n", xlab = "", ylab = "", xaxs = "i")
         abline(h=1:(length(tempVars)+length(catVars)), col = "grey")
         dumpVar3 = sapply(1:length(tempVars), function(i){
            pr = binVars[i]
            borders = diff(plotData[[tissue]][,pr])
            starts = which(borders > 0)
            ends = which(borders < 0)
            if(length(starts) > length(ends)){
               ends = max(pos)
            } else if(length(starts) < length(ends)){
               start = min(pos)
            }
            if(length(starts) > 0){
               segments(y0=i,y1 = i,
                        x0 = r2i(pos[starts]),
                        x1 = r2i(pos[ends]), lwd = 4, lend=1)
            }
         })
         dumpVar3 = sapply(1:length(catVars), function(i){
            points(r2i(pos), rep(i+length(binVars),length(pos)),
                   col = baseCols[plotData[[tissue]][,catVars[i]]], 
                   pch = "|", cex = 0.7)
         })
         axis(2,at = 1:(length(binVars)+length(catVars)), labels=c(binVars,catVars),
              tick=F, las = 1,
              mgp = c(0,0.5,0))
         rect(ybottom=par("usr")[3], ytop=par("usr")[4],
              xleft=r2i(holes[,1]),
              xright=r2i(holes[,2]), col="grey", border=NA)
         
         # numerical predictors
         tempVars = sapply(numVarsGrouped, function(y){
            y[y%in% colnames(plotData[[tissue]])]
         })
         tempVars = tempVars[sapply(tempVars,length)>0]
         tempCols = predicCols[names(tempVars)]
         plot(NA,xlim = range(r2i(pos)),
              ylim = c(1,length(tempVars)+1),xaxt ="n", yaxt = "n",
              xlab = "", ylab = "", xaxs = "i", yaxs = "i")
         text(x = min(r2i(pos)), y = 1:length(tempVars)+0.5, 
              labels = paste0(rev(names(tempVars))," "), 
              xpd = NA, adj = c(1, 0.5), col = darkCol(rev(tempCols), 2)[2,])
         dumpVar3 = sapply(1:length(tempVars), function(i){
            cl = darkCol(tempCols[i], length(tempVars[[i]]))
            for(j in length(tempVars[[i]]):1){
               temp = plotData[[tissue]][,tempVars[[i]][j]]
               temp = temp*0.85+0.05
               x = r2i(c(pos[1], pos, pos[length(pos)]))
               y = c(length(tempVars)-i+1,temp+length(tempVars)-i+1,length(tempVars)-i+1)
               polygon(x,y, border=NA, col = lightCol(cl[j], 100))
               lines(r2i(pos), temp+length(tempVars)-i+1, col = cl[j], lty = j,lwd = 2)
            }
         })
         
         rect(ybottom=par("usr")[3], ytop=par("usr")[4],
              xleft=r2i(holes[,1]),
              xright=r2i(holes[,2]), col="grey", border=NA)
         segments(y0 = 1:(length(tempVars)+1),col = "grey",
                  x0=r2i(max(pos)), x1 = r2i(max(pos))*1.1, xpd = NA)
         abline(h=1:length(tempVars), xlim = r2i(range(pos)))
         dumpVar3 = sapply(1:length(tempVars), function(i){
            if(length(tempVars[[i]])>1){
               labels = sapply(strsplit(tempVars[[i]], split="_"), 
                               function(x){paste(x[-1], collapse="_")})
               if(names(tempVars)[i] == "expression"){
                  labels = sapply(strsplit(tempVars[[i]], split="_"),
                                  function(x){x[1]})
               }
               if(names(tempVars)[i] %in% c("TF", "PolR2A")){
                  labels = tempVars[[i]]
               }
               labels = sapply(1:length(labels), function(k){
                  paste(c(rep("\n",k-1)," ",labels[k]),  sep="", collapse = "")
               })
               text(x=r2i(max(pos)), y=length(tempVars)-i+1.95,
                    labels=labels,cex = 0.6,font = 2,
                    xpd = NA, col=darkCol(tempCols[i], length(tempVars[[i]])),
                    adj=c(0,1))
            }
         })
         box()
         
         # add x-axis labels
         axis(1,at = r2i(seq(st, end, by=by)),
              labels=seq(st, end, by=by)/by,
              mgp = c(2,0.8,0))
         mtext(paste0("genomic position *",r[["label"]]), 
               side = 1, line=2, cex = 0.8)
         title(main = paste0(tissue, ", ",r["Gene"]), outer=T, 
               cex.main = 1.5, line = 1.5)
      }
      dev.off()
   }); cat('\n')
   
   
   # iterate through predictors and compare tissues
   print("predictor figures")
   dumpVar2 = sapply(names(numVarsTissue), function(pred){
      cat(pred, ' ')
      p = numVarsTissue[[pred]]
      pbyT = sapply(tissues, function(tissue){
         if(! any(p %in% colnames(plotData[[tissue]]))){
            return(NULL)
         } else{
            plotData[[tissue]][,p[p %in% colnames(plotData[[tissue]])]]
         }
      }, simplify = F)
      
      # plot tissues underneath each other
      png(paste0("fig/exomeVisualization/", r[["Gene"]], "_", pred, ".png"),
          height=2000, width=1400, res = 200)
      {
         par(mar = c(0,6.2,4,0), oma = c(4,0.5,0,0.5))
         layout(rbind(1,2,2,2,2,2,2,2,2,2,2,2))
         # gene annotation
         plot(NA, xlim = range(r2i(pos)), ylim = c(0,2),bty = "n",
              yaxt = "n",xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
         if(!is.na(r[["transcriptID"]])){
            abline(h=1)
            segments(x0=r2i(annot[["transcript"]]$V4), x1 = r2i(annot[["transcript"]]$V5), 
                     y0=1, lwd = 6, col = "black", lend=1)
            segments(x0=r2i(annot[["exon"]]$V4), x1 = r2i(annot[["exon"]]$V5), 
                     y0=1, lwd = 10, col = "black", lend=1)
            segments(x0=r2i(annot[["UTR"]]$V4), x1 = r2i(annot[["UTR"]]$V5), 
                     y0=1, lwd = 10, col = "blue", lend=1)
            segments(x0=r2i(annot[["CDS"]]$V4), x1 = r2i(annot[["CDS"]]$V5), 
                     y0=1, lwd = 10, col = "red", lend=1)
            if("stop_codon" %in% names(annot)){
               points(x=r2i(annot[["stop_codon"]]$V4), 
                      y=rep(1, nrow(annot[["stop_codon"]])), 
                      pch = "|", col = "red", cex = 2, xpd = NA)
            }
            if("start_codon" %in% names(annot)){
               points(x=r2i(annot[["start_codon"]]$V4), xpd = NA,
                      y=rep(1, nrow(annot[["start_codon"]])), 
                      pch = "|", col = "blue", cex = 1.5)
            }
            legend(x = r2i(max(pos)), y = 1, xpd = NA,ncol = 4,
                   legend = c("UTR", "exon", "start", "stop"), xjust = 1,yjust = 0,
                   y.intersp=0.9,bty = "n",cex = 1.5,
                   pch = c(NA, NA, "|", "|"), col = c(NA, NA, "blue", "red"),
                   fill = c("blue", "red", NA, NA), border=NA)
            legend(x = r2i(min(pos)), y = 1, xpd = NA,
                   xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",cex = 1.5,
                   legend = paste(ifelse(r[["strand"]]=="+", yes="---->", no = "<----"),
                                  "reading direction"),border=NA)
         }
         # predictors
         plot(NA,xlim = range(r2i(pos)), 
              ylim = c(1,length(tissues)+1),xaxt ="n", yaxt = "n",
              xlab = "", ylab = "", xaxs = "i", yaxs = "i")
         text(x = min(r2i(pos)), y = 1:length(tissues)+0.5, 
              labels = paste0(rev(tissues)," "), cex = 1.5,
              xpd = NA, adj = c(1, 0.5), col = darkCol(rev(tissueCols), 2)[2,])
         tempCol = darkCol(tissueCols, length(p))
         dumpVar3 = sapply(1:length(tissues), function(i){
            for(j in 1:length(p)){
               if(is.null(pbyT[[tissues[i]]]) | 
                  !p[j] %in% colnames(pbyT[[tissues[i]]])){next}
               temp = pbyT[[tissues[i]]][,p[j]]
               # rang = pRanges[[p[j]]]
               # temp = (temp-rang[1])/(rang[2]-rang[1])
               # temp = temp*0.90+0.05
               temp = temp*0.65+0.05
               x = r2i(c(pos[1], pos, pos[length(pos)]))
               y = c(length(tissues)-i+1,temp+length(tissues)-i+1,length(tissues)-i+1)
               lines(r2i(pos), lwd = 3,
                     temp+length(tissues)-i+1, col = tempCol[j,tissues[i]],
                     lty = j)
               polygon(x,y, border=NA, col = lightCol(tempCol[j,tissues[i]], 100))
            }
            muts = plotData[[tissues[i]]]$mutated
            abline(h=length(tissues)-i+1.85, col = "grey")
            if(any(muts)){
               # segments(x0 = r2i(pos[muts]),
               #          y0=length(tissues)-i+1, y1 = length(tissues)-i+2)
               segments(x0 = r2i(pos[muts]),
                        y0=length(tissues)-i+1.75, y1 = length(tissues)-i+1.95)
            }
         })
         rect(ybottom=par("usr")[3], ytop=par("usr")[4], 
              xleft=r2i(holes[,1]),
              xright=r2i(holes[,2]), col="grey", border=NA)
         axis(1,at = r2i(seq(st, end, by=by)), 
              labels=seq(st, end, by=by)/by,
              mgp = c(2,0.8,0), cex.axis = 1.5)
         mtext(paste0("genomic position *",r[["label"]]), side = 1, line=2)
         if(length(p)>1){
            legend(x = r2i(max(pos)), y = length(tissues)+1, xpd = NA,ncol = 2,
                   legend = p, xjust = 1,yjust = 0,y.intersp=0.9,bty = "n",cex = 1.5,
                   lty = 1:length(p), lwd = 3, col = darkCol("grey", length(p)))
         }
         legend(x = r2i(min(pos)), y = length(tissues)+1, xpd = NA,
                legend = "mutation", xjust = 0,yjust = 0,y.intersp=0.9,x.intersp=0.9,
                bty = "n",cex = 1.5,lty = 1, col = "grey",  seg.len = 1.5, merge = T)
         legend(x = r2i(min(pos)), y = length(tissues)+1, xpd = NA,
                legend = "", xjust = 0,yjust = 0,y.intersp=0.9,x.intersp=0.9,
                bty = "n",cex = 1.5,pch = "|")
         abline(h=1:length(tissues))
         title(main = paste0(pred, ", ",r["Gene"]), line = -2, 
               cex.main = 1.5, outer = T)
         box()
      }
      dev.off()
   }); cat('\n')
   
   
   
   # compare predictions across tissues
   # all predictors
   print("predictions")
   png(paste0("fig/exomeVisualization/", r[["Gene"]], "_predictions.png"),
       height=2000, width=1400, res = 200)
   {
      par(mar = c(0,6,0,0), oma = c(3.1,0.1,3,0.1))
      layout(cbind(c(1,rep(1:length(predVars)+1, each = 3))))
      plot(NA, xlim = range(r2i(pos)), ylim = c(-2,2),bty = "n",
           yaxt = "n",xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
      if(!is.na(r[["transcriptID"]])){
         abline(h=1)
         segments(x0=r2i(annot[["transcript"]]$V4), x1 = r2i(annot[["transcript"]]$V5), 
                  y0=1, lwd = 6, col = "black", lend=1)
         segments(x0=r2i(annot[["exon"]]$V4), x1 = r2i(annot[["exon"]]$V5), 
                  y0=1, lwd = 10, col = "black", lend=1)
         segments(x0=r2i(annot[["UTR"]]$V4), x1 = r2i(annot[["UTR"]]$V5), 
                  y0=1, lwd = 10, col = "blue", lend=1)
         segments(x0=r2i(annot[["CDS"]]$V4), x1 = r2i(annot[["CDS"]]$V5), 
                  y0=1, lwd = 10, col = "red", lend=1)
         if("stop_codon" %in% names(annot)){
            points(x=r2i(annot[["stop_codon"]]$V4), 
                   y=rep(1, nrow(annot[["stop_codon"]])), 
                   pch = "|", col = "red", cex = 2)
         }
         if("start_codon" %in% names(annot)){
            points(x=r2i(annot[["start_codon"]]$V4), 
                   y=rep(1, nrow(annot[["start_codon"]])), 
                   pch = "|", col = "blue", cex = 1.5)
         }
         legend(x = r2i(max(pos)), y = 1.2, xpd = NA,ncol = 4,
                legend = c("UTR", "exon", "start", "stop"), xjust = 1,yjust = 0,
                y.intersp=0.9,bty = "n",cex = 1.2,
                pch = c(NA, NA, "|", "|"), col = c(NA, NA, "blue", "red"),
                fill = c("blue", "red", NA, NA), border=NA)
         legend(x = r2i(min(pos)), y = 1.2, xpd = NA,cex = 1.2,
                xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",
                legend = paste(ifelse(r[["strand"]]=="+", yes="---->", no = "<----"),
                               "reading direction"),border=NA)
      }
      legend(x = r2i(max(pos)), y = -2.5, xpd = NA,cex = 1.2,
             xjust = 1,yjust = 0,y.intersp=0.9,bty = "n",border=NA,
             legend = "mutation", pch = "|")
      
      # predictors
      dumpVar2 = sapply(predVars, function(p){
         # pRange = range(sapply(plotData, function(x){range(x[,p])}))
         plot(NA,xlim = range(r2i(pos)),
              ylim = c(1,length(tissues)+1),xaxt ="n", yaxt = "n",
              xlab = "", ylab = "", xaxs = "i", yaxs = "i")
         text(x = min(r2i(pos)), y = 1:length(tissues)+0.5, 
              labels = paste0(rev(tissues)," "), 
              xpd = NA, adj = c(1, 0.5), col = rev(tissueCols))
         title(ylab = p, line = 5, cex.lab = 1.2)
         for(i in 1:length(tissues)){
            temp = plotData[[tissues[i]]][,p]
            # temp = (temp-pRange[1])/(pRange[2]-pRange[1])
            # temp = (temp-min(temp))/(max(temp)-min(temp))
            temp = temp*0.95+0.025
            x = r2i(c(pos[1], pos, pos[length(pos)]))
            y = c(length(tissues)-i+1,temp+length(tissues)-i+1,length(tissues)-i+1)
            polygon(x,y, col = tissueCols[tissues[i]], border=NA)
            muts = plotData[[tissues[i]]]$mutated
            if(any(muts)){
               segments(x0 = r2i(pos[muts]), 
                        y0=length(tissues)-i+2, y1 = length(tissues)-i+1)
            }
         }
         rect(ybottom=par("usr")[3], ytop=par("usr")[4], 
              xleft=r2i(holes[,1]),
              xright=r2i(holes[,2]), col="grey", border=NA)
         abline(h = 1:length(tissues))
         box()
         abline(h = c(1,length(tissues)+1), lwd = 1, 
                xlim = c(0,max(r2i(pos))), xpd = NA)
      })
      axis(1,at = r2i(seq(st, end, by=by)), 
           labels=seq(st, end, by=by)/by,
           mgp = c(2,0.8,0))
      mtext(paste0("genomic position *",r[["label"]]), side = 1, line=2)
      title(main = r["Gene"],  outer = T, line = 1.5)
   }
   dev.off()
   
   # only three predictors
   print("predictions")
   png(paste0("fig/exomeVisualization/", r[["Gene"]], "_predictions2.png"),
       height=2000, width=1400, res = 200)
   {
      par(mar = c(0,6,0,0), oma = c(3.1,0.1,3,0.1))
      tempVars = c("RFpreds", "contextPreds", "combOdds")
      layout(cbind(c(1,rep(1:length(tempVars)+1, each = 3))))
      plot(NA, xlim = range(r2i(pos)), ylim = c(-2,2),bty = "n",
           yaxt = "n",xaxt = "n", ylab = "", xaxs = "i", yaxs = "i")
      if(!is.na(r[["transcriptID"]])){
         abline(h=1)
         segments(x0=r2i(annot[["transcript"]]$V4), x1 = r2i(annot[["transcript"]]$V5), 
                  y0=1, lwd = 6, col = "black", lend=1)
         segments(x0=r2i(annot[["exon"]]$V4), x1 = r2i(annot[["exon"]]$V5), 
                  y0=1, lwd = 10, col = "black", lend=1)
         segments(x0=r2i(annot[["UTR"]]$V4), x1 = r2i(annot[["UTR"]]$V5), 
                  y0=1, lwd = 10, col = "blue", lend=1)
         segments(x0=r2i(annot[["CDS"]]$V4), x1 = r2i(annot[["CDS"]]$V5), 
                  y0=1, lwd = 10, col = "red", lend=1)
         if("stop_codon" %in% names(annot)){
            points(x=r2i(annot[["stop_codon"]]$V4), 
                   y=rep(1, nrow(annot[["stop_codon"]])), 
                   pch = "|", col = "red", cex = 2)
         }
         if("start_codon" %in% names(annot)){
            points(x=r2i(annot[["start_codon"]]$V4), 
                   y=rep(1, nrow(annot[["start_codon"]])), 
                   pch = "|", col = "blue", cex = 1.5)
         }
         legend(x = r2i(max(pos)), y = 1.2, xpd = NA,ncol = 4,
                legend = c("UTR", "exon", "start", "stop"), xjust = 1,yjust = 0,
                y.intersp=0.9,bty = "n",cex = 1.2,
                pch = c(NA, NA, "|", "|"), col = c(NA, NA, "blue", "red"),
                fill = c("blue", "red", NA, NA), border=NA)
         legend(x = r2i(min(pos)), y = 1.2, xpd = NA,cex = 1.2,
                xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",
                legend = paste(ifelse(r[["strand"]]=="+", yes="---->", no = "<----"),
                               "reading direction"),border=NA)
      }
      legend(x = r2i(min(pos)), y = -1.5, xpd = NA,cex = 1.2,
             xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",border=NA,
             legend = "mutation", lty = 1, col = "grey", seg.len = 1.5, merge = T)
      legend(x = r2i(min(pos)), y = -1.5, xpd = NA,cex = 1.2,
             xjust = 0,yjust = 0,y.intersp=0.9,bty = "n",border=NA,
             legend = "", pch = "|")
      # predictors
      dumpVar2 = sapply(tempVars, function(p){
         pRange = range(sapply(plotData, function(x){range(x[,p])}))
         plot(NA,xlim = range(r2i(pos)),
              ylim = c(1,length(tissues)+1),xaxt ="n", yaxt = "n",
              xlab = "", ylab = "", xaxs = "i", yaxs = "i")
         text(x = min(r2i(pos)), y = 1:length(tissues)+0.5, 
              labels = paste0(rev(tissues)," "), 
              xpd = NA, adj = c(1, 0.5), col = rev(darkCol(tissueCols, 2)[2,]))
         title(ylab = p, line = 5, cex.lab = 1.2)
         for(i in 1:length(tissues)){
            temp = plotData[[tissues[i]]][,p]
            # if(p == "RFpreds"){
            #    temp = (temp-pRange[1])/(pRange[2]-pRange[1])
            # } else{
            #    temp = (temp-min(temp))/(max(temp)-min(temp))
            # }
            temp = temp*0.7+0.025
            x = r2i(c(pos[1], pos, pos[length(pos)]))
            y = c(length(tissues)-i+1,temp+length(tissues)-i+1,length(tissues)-i+1)
            polygon(x,y, col = darkCol(tissueCols[tissues[i]], 2)[2,], border=NA)
            muts = plotData[[tissues[i]]]$mutated
            abline(h=length(tissues)-i+1.75, col = "grey")
            if(any(muts)){
               segments(x0 = r2i(pos[muts]), 
                        y0=length(tissues)-i+1.6, y1 = length(tissues)-i+1.9)
            }
         }
         rect(ybottom=par("usr")[3], ytop=par("usr")[4], 
              xleft=r2i(holes[,1]),
              xright=r2i(holes[,2]), col="grey", border=NA)
         abline(h = 1:length(tissues))
         box()
         abline(h = c(1,length(tissues)+1), lwd = 1, 
                xlim = c(0,max(r2i(pos))), xpd = NA)
      })
      axis(1,at = r2i(seq(st, end, by=by)), 
           labels=seq(st, end, by=by)/by,
           mgp = c(2,0.8,0))
      mtext(paste0("genomic position *",r[["label"]]), side = 1, line=2)
      title(main = r["Gene"],  outer = T, line = 1.5)
   }
   dev.off()
})
#####

print("done")


