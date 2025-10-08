#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmalo
#SBATCH --output=output/prepareHealthyMuts_o_%j
#SBATCH --error=output/prepareHealthyMuts_e_%j
#SBATCH -w beyer-n04

module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/02_prepareMutations/03_healthyTissues/01_healthyTissuesWEXandWGS.R' 
