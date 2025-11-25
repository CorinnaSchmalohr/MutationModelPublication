#!/bin/bash -l
#SBATCH --account=cschmal1
#SBATCH --output=output/mapHealthyMutsWEX_o_%j
#SBATCH --error=output/mapHealthyMutsWEX_e_%j
#SBATCH --array=1-29

module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/05_healthyTissues_WEX.R'  ${SLURM_ARRAY_TASK_ID}
