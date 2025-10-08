#!/bin/bash -l
#SBATCH --cpus-per-task=56
#SBATCH --mem 250G
#SBATCH --account=cschmalo
#SBATCH --output=output/09b_generalModelEvaluation_WGS_o_%j
#SBATCH --error=output/09b_generalModelEvaluation_WGS_e_%j
#SBATCH -w beyer-n04

module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose "scripts/05_analysis/02_WGSdata/02_generalModelEvaluation_WGS.R"
