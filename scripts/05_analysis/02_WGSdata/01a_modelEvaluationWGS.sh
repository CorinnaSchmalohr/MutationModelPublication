#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmal1
#SBATCH --output=output/01c_modelEvaluationWGS_e_%j
#SBATCH --error=output/01c_modelEvaluationWGS_o_%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/02_WGSdata/01a_modelEvaluationWGS.R' 

