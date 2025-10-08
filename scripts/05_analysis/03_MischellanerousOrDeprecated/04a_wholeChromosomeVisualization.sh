#!/bin/bash -l
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/WEvis_o-%j
#SBATCH --error=output/WEvis_e-%j

source /data/public/apapada1/Conda/anaconda/bin/activate
conda activate r4-base
Rscript --vanilla --verbose 'scripts/05_analysis/04a_wholeChromosomeVisualization.R' 
