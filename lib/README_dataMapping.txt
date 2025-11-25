The predictors were mapped to query genomic positions (e.g., TPs and TNs) using a streamlined approach. 
For each predictor and each tissue, we created one file encompassing genome-wide values for this predictor, 
based on pre-processed and collated data from multiple files for each tissue. Depending on the predictor, 
we used either GR/GRL objects or bigWig files, with exception of GC content and the effect score, 
which are computed for each position in-place. To create the dataset for model training and/or prediction, 
we simply read out the value corresponding to each query position (i.e., TP or TN positions) 
for each tissue from each GR or bigWig, using a custom R function. This function takes as input 
the file location of the bed file with the query positions, as well as a table indicating which 
predictors should be mapped in which way, thus acting as a configuration file. This table includes
the file locations of the pre-processed predictor to be used, as well as further options for mapping.
Relevant here are the following parameters: "range" indicates the size of the window 
around the query position that should be considered when reading out a predictor, 
from 1bp up to 1Mbp, while an "N/A" indicates only a 1bp window. The option "measure" 
defines the manner with which a predictor should be read out. Possible options are 
"distance" (distance from query position to closest feature in the GR object in bps, 
used for example to compute distance to telomeres and centromeres), 
"ifany" (indicating whether the query position overlaps any feature in the GR at all,
resulting in a binary 0/1 predictor), "mean" (mean predictor value across a range 
around the positions, with regions that are not covered by the reference predictor file 
being ignored), "mean0" (same as "mean", only that non-covered regions are counted as 0 
towards the mean), and "nHits" (number of GR features within the window around the query position). 
Finally, the option "transform" indicates whether the predictor should be further 
transformed by either taking the square root ("sqrt"), the logarithm ("log"), or used as-is (NA). 
The transform options were chosen to approximate the predictor distributions to a gaussian 
distribution and make them comparable across tissues and file sources.
The mapPredictors function extracts values from BigWigs using BigWigAverageOverBed, 
the "range" option is simply passed on to -sampleAroundCenter. GR objects were read out 
by creating a GR object from the query positions as well and extending the ranges in either 
direction according to the by range option. For the output options "mean" or "mean0", we used 
the function findOverlaps to find matches between query positions and predictor GR and took the 
mean of all overlap hits, counting non-covered positions as either 0 or NA for "mean0" or "mean", 
respectively. When taking the mean, overlap width is taken into account (i.e., if only part of the 
range around a variant was covered in the GR object, only the regions which were actually overlapping 
were taken into account, the rest was considered non-covered). For output "measure" options "ifany", 
"nHits", or "distance", we used the output of the functions findOverlaps, countOverlaps, or distanceToNearest, respectively.

