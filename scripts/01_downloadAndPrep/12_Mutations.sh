# TCGA Pancan  mutations
wget -O data/rawdata/pancan/mc3.v0.2.8.PUBLIC.maf.gz \
  http://api.gdc.cancer.gov/data/1c8cfe5f-e52d-41ba-94da-f15ea1337efc


# ICGC PCAWG
wget -O data/rawdata/pancan/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz \
https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz


# SomaMutDB V1.4 (January 2024)
mkdir data/rawdata/SomaMutDB_hg19/
wget --no-check-certificate -P data/rawdata/SomaMutDB_hg19/ \
https://vijglab.einsteinmed.edu/static/vcf/download.txt \
https://vijglab.einsteinmed.edu/static/vcf/hg19_34.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_101.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_24.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_3.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_45.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_45.tar.gz
https://vijglab.einsteinmed.edu/static/vcf/hg19_10.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_43.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_51.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_42.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_50.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_52.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_2.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_46.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_48.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_54.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_16.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_30.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_8.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_28.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_9.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_37.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_27.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_6.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_11.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_4.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_38.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_53.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_106.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_104.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_105.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_103.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_102.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_18.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_41.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_17.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_12.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_33.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_20.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_32.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_47.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_22.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_35.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_19.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_1.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_23.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_7.tar.gz \
https://vijglab.einsteinmed.edu/static/vcf/hg19_13.tar.gz
cd data/rawdata/SomaMutDB_hg19
flist=$(ls *.tar.gz)
for f in ${flist}; do
echo $f
tar -xzf $f
done
