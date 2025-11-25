#!/bin/bash -l
#SBATCH --account=cschmal1
#SBATCH --output=output/mapHealthyMutsWGS_o_%j
#SBATCH --error=output/mapHealthyMutsWGS_e_%j
#SBATCH --array=1-28

module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/03_mapPredictors/05_healthyTissues_WGS.R'  ${SLURM_ARRAY_TASK_ID}
