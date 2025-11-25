#!/bin/bash -l
#SBATCH --job-name=genome_Mapping
#SBATCH --cpus-per-task=32
#SBATCH --account=ypaul1
#SBATCH --output=./tmp/output/genomeRF_o-%A_%a.txt
#SBATCH --error=./tmp/output/genomeRF_e-%A_%a.txt
#SBATCH --array=18,17
#SBATCH --mem=150
#SBATCH --time=100:00:00



# Initialize conda
source ~/.bashrc  # Or your conda init file
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R

CHR_IDX=$SLURM_ARRAY_TASK_ID


# Step 2: Process parts 1-20
echo "Processing parts 1-20 for chromosome $CHR_IDX"
Rscript --vanilla '/cellfile/cellnet/MutationModel/scripts/03_mapPredictors/02_WGSdata/2025_11_17_04e_New_PartsThatRemain.R' $CHR_IDX "1-20"

echo "Completed processing for chromosome $CHR_IDX"
