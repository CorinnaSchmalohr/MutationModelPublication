#!/bin/bash -l
#SBATCH --cpus-per-task=32
#SBATCH --account=cschmal1
#SBATCH --output=output/crossTmuttype_o-%j
#SBATCH --error=output/crossTmuttype_e-%j


module unload R-3.5.1
module load R-4.1.2
Rscript --vanilla --verbose 'scripts/05_analysis/01_WESdata/04b_crossTissueApplicationMuttype.R' 
