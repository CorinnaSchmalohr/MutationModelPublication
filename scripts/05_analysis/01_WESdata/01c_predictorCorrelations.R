library(readxl)
source("scripts/05_analysis/00_NamesAndColors.R")
plotEnding =  "_20240925"

# prepare predictor annotation #####
tab = read_xlsx("data/rawdata/dataMappingAlltissues.xlsx", 
sheet="allTissues", col_names=T)
ranges = c("1Mb" = 500000, "100kb" = 50000,"10kb" = 5000, 
           "1kb" = 500, "100bp" = 50, "10bp" = 5, "1bp" = 0)
tab$NA. = NULL
tab[tab == "NA"] = NA
# for predictors where we want multiple ranges, expand table 
tab = apply(tab,1,function(x){
  if(is.na(x["range"])){
    return(x)
  } else{
    rangeWindow = strsplit(x["range"],";")[[1]]
    rangeInds = which(names(ranges) == rangeWindow[2]):which(names(ranges) == rangeWindow[1])
    subRanges = ranges[rangeInds]
    t(sapply(names(subRanges), function(r,y){
      y["range"] = subRanges[r]
      if(length(subRanges)>1){
        y["Name"] = paste0(y["Name"]," ",r)
        y["abbreviation"] = paste0(y["abbreviation"],"_",r)
      }
      return(y)
    },y=x))
  }
})
tab = do.call(rbind, tab)
tab = as.data.frame(tab)
predictors = tab$abbreviation[tab$tissueSpecific == "yes"]
######

# compute predictor correlations #####
predictorCors = sapply(predictors, function(predictor){
  print(predictor)
  # for each tissue, collect the predictor info from its own tissue 
  # as well as the cross-tissue mappings in the other tissues
  predictorCollection = sapply(tissues, function(predtissue){
    # print(predtissue)
    if(is.na((tab[tab$abbreviation == predictor,predtissue]))){
      return(NA) # in case predictor was not available for this tissue
    }
    predictorTissue = sapply(tissues, function(muttissue){
      # print(muttissue)
      if(muttissue == predtissue){
        load(paste0("data/MutTables/exomeTrainData/", 
                    muttissue, "_Muts_mapped_processed.RData"))
      } else{
        load(paste0("data/MutTables/exomeTrainDataCrossTissue/", 
                    muttissue, "With",predtissue,"Preds_processed", ".RData")) # dat, datchroms
      }
      return(dat[,predictor])
    }, simplify = F)
    return(unlist(predictorTissue))
  })
  temp = as.data.frame(predictorCollection)
  res = cor(temp, use = "pair")
  return(res)
}, simplify = F)
save(predictorCors, file = "data/MutTables/exomeTrainDataCrossTissue/predictorCors.R")
######

# visualize #####
# load("data/MutTables/exomeTrainDataCrossTissue/predictorCors.R")

library(ggplot2)
library(colorspace)
plotDat = reshape2::melt(predictorCors)
plotDat$rounded = format(plotDat$value, digits = 2)
plotDat$rounded[is.na(plotDat$value)] = ""
plotDat$L1 = stringr::str_wrap(p2P[plotDat$L1], width = 12)
plotDat$L1 = factor(plotDat$L1, levels = stringr::str_wrap(p2P, width = 12))

ggplot(plotDat, aes(x=Var1, y = Var2, fill = value))+
  geom_tile() +
  # scale_fill_gradient2(low = "#075AFF",
  #                      mid = "#FFFFCC",
  #                      high = "#FF0000") +
  scale_fill_continuous_divergingx(palette = 'RdBu') + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y=element_text(size = 5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=5, lineheight=0.8, margin = margin(0,0,0,0)),
        panel.spacing = unit(0.1, "lines"))+
  labs(y = "", x = "") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  # geom_text(aes(label = rounded), color = "white", size = 2.5)+
  coord_fixed() +
  guides(fill = guide_colourbar(title = "Pearson R",
                                barwidth = 0.5,
                                barheight = 10))+
  facet_wrap(~L1, nrow = 6) 
ggsave(paste0("fig/predictorCrossTissueCorrelations", plotEnding, ".png"),
       width = 30, height = 15,units = "cm")


#####


# compute predictor spearman correlations #####
predictorCorsSpearman = sapply(predictors, function(predictor){
  print(predictor)
  # for each tissue, collect the predictor info from its own tissue 
  # as well as the cross-tissue mappings in the other tissues
  predictorCollection = sapply(tissues, function(predtissue){
    # print(predtissue)
    if(is.na((tab[tab$abbreviation == predictor,predtissue]))){
      return(NA) # in case predictor was not available for this tissue
    }
    predictorTissue = sapply(tissues, function(muttissue){
      # print(muttissue)
      if(muttissue == predtissue){
        load(paste0("data/MutTables/exomeTrainData/", 
                    muttissue, "_Muts_mapped_processed.RData"))
      } else{
        load(paste0("data/MutTables/exomeTrainDataCrossTissue/", 
                    muttissue, "With",predtissue,"Preds_processed", ".RData")) # dat, datchroms
      }
      return(dat[,predictor])
    }, simplify = F)
    return(unlist(predictorTissue))
  })
  temp = as.data.frame(predictorCollection)
  res = cor(temp, use = "pair", method = "spearman")
  return(res)
}, simplify = F)
save(predictorCorsSpearman, file = "data/MutTables/exomeTrainDataCrossTissue/predictorCorsSpearman.R")
######


# visualize #####
# load("data/MutTables/exomeTrainDataCrossTissue/predictorCorsSpearman.R")

library(ggplot2)
library(colorspace)
plotDat = reshape2::melt(predictorCorsSpearman)
plotDat$rounded = format(plotDat$value, digits = 2)
plotDat$rounded[is.na(plotDat$value)] = ""
plotDat$L1 = stringr::str_wrap(p2P[plotDat$L1], width = 12)
plotDat$L1 = factor(plotDat$L1, levels = stringr::str_wrap(p2P, width = 12))

ggplot(plotDat, aes(x=Var1, y = Var2, fill = value))+
  geom_tile() +
  # scale_fill_gradient2(low = "#075AFF",
  #                      mid = "#FFFFCC",
  #                      high = "#FF0000") +
  scale_fill_continuous_divergingx(palette = 'RdBu') + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
        axis.text.y=element_text(size = 5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_text(size=5, lineheight=0.8, margin = margin(0,0,0,0)),
        panel.spacing = unit(0.1, "lines"))+
  labs(y = "", x = "") +  
  scale_x_discrete(labels=t2T) + 
  scale_y_discrete(labels = t2T)+
  # geom_text(aes(label = rounded), color = "white", size = 2.5)+
  coord_fixed() +
  guides(fill = guide_colourbar(title = "Pearson R",
                                barwidth = 0.5,
                                barheight = 10))+
  facet_wrap(~L1, nrow = 6) 
ggsave(paste0("fig/predictorCrossTissueSpearmanCorrelations", plotEnding, ".png"),
       width = 30, height = 15,units = "cm")


#####