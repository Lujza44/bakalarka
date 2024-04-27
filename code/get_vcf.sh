#!/bin/bash

ORIGINAL_DIR=$(pwd)
cd data/input

wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/00-common_all.vcf.gz
gunzip 00-common_all.vcf.gz
sudo apt-get install tabix
bgzip -c 00-common_all.vcf > 00-common_all.vcf.gz
bcftools index 00-common_all.vcf.gz

cd "$ORIGINAL_DIR"