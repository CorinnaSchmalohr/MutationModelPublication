#!/bin/bash -l
#SBATCH --cpus-per-task=56
#SBATCH --mem 250G
#SBATCH --account=cschmalo
#SBATCH --output=output/RFWGS_o-%j
#SBATCH --error=output/RFWGS_e-%j
#SBATCH --array=8
#SBATCH --partition=all

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/02_WGSdata/01_RF_WGSdata.R' ${SLURM_ARRAY_TASK_ID}
