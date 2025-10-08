#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --account=cschmal1
#SBATCH --output=output/codExons_noUTR_filtered_o-%j
#SBATCH --error=output/codExons_noUTR_filtered_e-%j

# translate gtf into bed
cat data/processedData/codExons_noUTR_filtered.gtf |
  awk 'BEGIN{OFS="\t"}; {$4=$4-3; $5=$5+2; print $1,$4,$5, $1" "$4" "$5" "$8}' | 
  bedtools getfasta -fi data/rawdata/genome/GRCh37.primary_assembly.genome.fa \
-bed - \
-name -tab > data/processedData/codExons_noUTR_filtered.withsequences.bed

echo 'done'
