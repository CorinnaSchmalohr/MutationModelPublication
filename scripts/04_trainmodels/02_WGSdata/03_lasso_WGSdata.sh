#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --output=output/lassoWGS_o-%j
#SBATCH --error=output/lassoWGS_e-%j
#SBATCH --array=1-8


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/02_WGSdata/03_lasso_WGSdata.R' ${SLURM_ARRAY_TASK_ID}
