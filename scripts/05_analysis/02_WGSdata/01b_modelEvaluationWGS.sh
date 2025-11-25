#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --mem=80G  
#SBATCH --account=cschmal1
#SBATCH --output=output/WGSanalysisPlots_o-%j
#SBATCH --error=output/WGSanalysisPlots_e-%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/02_WGSdata/01b_modelEvaluationWGS.R' 
