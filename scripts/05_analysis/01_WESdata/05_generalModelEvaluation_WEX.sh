#!/bin/bash -l
#SBATCH --cpus-per-task=14
#SBATCH --account=cschmalo
#SBATCH --output=output/generalModelEvaluation_WEX_o_%j
#SBATCH --error=output/generalModelEvaluation_WEX_e_%j

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose "scripts/05_analysis/01_WESdata/05_generalModelEvaluation_WEX.R"
