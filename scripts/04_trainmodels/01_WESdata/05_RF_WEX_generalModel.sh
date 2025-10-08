#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=cschmalo
#SBATCH --output=output/05_RF_WEX_generalModel_o_%j
#SBATCH --error=output/05_RF_WEX_generalModel_e_%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/01_WESdata/05_RF_WEX_generalModel.R' 
