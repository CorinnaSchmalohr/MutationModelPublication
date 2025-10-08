#!/bin/bash -l
#SBATCH --cpus-per-task=56
#SBATCH --account=cschmal1
#SBATCH --output=output/glmWGS_o-%j
#SBATCH --error=output/glmWGS_e-%j
#SBATCH --partition=all
#SBATCH --array=1-8

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/02_WGSdata/02_glm_WGSdata.R' ${SLURM_ARRAY_TASK_ID}
