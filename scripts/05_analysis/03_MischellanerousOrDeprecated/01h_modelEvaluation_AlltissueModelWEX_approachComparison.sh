#!/bin/bash -l
#SBATCH --cpus-per-task=14
#SBATCH --account=cschmal1
#SBATCH --output=output/01h_modelEvaluation_AlltissueModelWEX_o_%j
#SBATCH --error=output/01h_modelEvaluation_AlltissueModelWEX_e_%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/01h_modelEvaluation_AlltissueModelWEX_approachComparison.R' 
