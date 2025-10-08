#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --output=output/glm_o-%j
#SBATCH --error=output/glm_e-%j
#SBATCH --array=1-10

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/01_WESdata/02_glm.R' ${SLURM_ARRAY_TASK_ID}
