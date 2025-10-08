#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=nbundsc1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=nbundsc1@uni-koeln.de
#SBATCH --output=output/glm_SelfSelection_o-%j
#SBATCH --error=output/glm_SelfSelection_e-%j
#SBATCH --array=1-10


module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/10_glm_FeatureSource_SelfSelection.R' ${SLURM_ARRAY_TASK_ID}
