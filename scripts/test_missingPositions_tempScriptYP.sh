#!/bin/bash -l
#SBATCH --job-name=test
#SBATCH --cpus-per-task=56
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/RF_o-%A_%a.txt
#SBATCH --error=./tmp/output/RF_e-%A_%a.txt
#SBATCH --mem=150
#SBATCH --time=100:00:00


# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R


# Run R script
Rscript --vanilla '/cellfile/cellnet/MutationModel/scripts/test_missingPositions_tempScriptYP.R'

