#!/bin/bash

set -e

ORIGINAL_DIR=$(pwd)
cd data/input

VCF_FILE="00-common_all.vcf"
VCF_GZ_FILE="$VCF_FILE.gz"

if [ -f "$VCF_FILE" ] || [ -f "$VCF_GZ_FILE" ]; then
    echo "VCF file already exists. Skipping download."
else
    wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/$VCF_GZ_FILE
    gunzip $VCF_GZ_FILE
    sudo apt-get install -y tabix
    bgzip -c $VCF_FILE > $VCF_GZ_FILE
    bcftools index $VCF_GZ_FILE
fi

cd "$ORIGINAL_DIR"
