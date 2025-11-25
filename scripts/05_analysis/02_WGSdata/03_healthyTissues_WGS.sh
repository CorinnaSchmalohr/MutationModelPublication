#!/bin/bash -l
#SBATCH --cpus-per-task=56
#SBATCH --mem 250G
#SBATCH --account=cschmalo
#SBATCH --output=output/healthyTissues_WGS_e_%j
#SBATCH --error=output/healthyTissues_WGS_o_%j
#SBATCH -w beyer-n04


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/02_WGSdata/03_healthyTissues_WGS.R' 

