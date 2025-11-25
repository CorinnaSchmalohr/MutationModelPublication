#!/bin/bash -l
#SBATCH --cpus-per-task=12
#SBATCH --account=cschmal1
#SBATCH --mail-type=END,FAIL,REQUEUE
#SBATCH --mail-use=cschmal1@uni-koeln.de
#SBATCH --output=output/prepareHealthyMuts_o_%j
#SBATCH --error=output/prepareHealthyMuts_e_%j
#SBATCH --array=1-26
#
module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/02_prepareMutations/05b_healthyTissues_WGS.R'  ${SLURM_ARRAY_TASK_ID}
