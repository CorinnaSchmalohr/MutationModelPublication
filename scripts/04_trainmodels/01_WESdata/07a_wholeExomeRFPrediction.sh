#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=cschmalo
#SBATCH --output=output/exomeRF_o-%j
#SBATCH --error=output/exomeRF_e-%j
#SBATCH --array=1-10

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/01_WESdata/07a_wholeExomeRFPrediction.R' ${SLURM_ARRAY_TASK_ID}
