#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmal1
#SBATCH --output=output/prepareWGSmuts_o-%j
#SBATCH --error=output/prepareWGSmuts_e-%j
#SBATCH --array=8

module unload R-3.5.1 
module load R-4.1.2
# source /data/public/cschmalo/anaconda3/bin/activate 
# conda activate test
Rscript --vanilla --verbose 'scripts/02_prepareMutations/02_WGSdata/02_prepareWGSPanCanMuts.R'  ${SLURM_ARRAY_TASK_ID}

