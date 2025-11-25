#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmal1
#SBATCH --output=output/WGSandWEX_o-%j
#SBATCH --error=output/WGSandWEX_e-%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/02_WGSdata/01c_WGSandWEXcrossApplication.R'
