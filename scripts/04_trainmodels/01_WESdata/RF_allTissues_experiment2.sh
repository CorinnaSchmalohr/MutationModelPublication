#!/bin/bash -l
#SBATCH --cpus-per-task=28
#SBATCH --account=cschmalo
#SBATCH --output=output/RFtissueCombination_o_%j
#SBATCH --error=output/RFtissueCombination_e_%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/04_trainmodels/09c_RF_allTissues_experiment2.R' 
