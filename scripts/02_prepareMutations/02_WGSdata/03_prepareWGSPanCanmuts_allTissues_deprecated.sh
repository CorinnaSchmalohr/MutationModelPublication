#!/bin/bash -l
#SBATCH --cpus-per-task=16
#SBATCH --account=cschmal1
#SBATCH --output=output/prepareAllTissueWGSmuts_o-%j
#SBATCH --error=output/prepareAllTissueWGSmuts_e-%j

module unload R-3.5.1 
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/02_prepareMutations/02_WGSdata/03_prepareWGSPanCanmuts_allTissues.R' 

