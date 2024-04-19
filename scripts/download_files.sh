#!/bin/bash

###### SETTINGS ######
# ENSEMBL release version
release=111
# set directory to store reference files
ref_dir=/path/to/reference_dir

###### START DOWNLOAD ######
# download reference files from ENSEMBL FTP
# below are the links for mouse reference files
# for other species, please browse the ENSEMBL FTP website
cd $ref_dir
wget https://ftp.ensembl.org/pub/release-${release}/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
wget https://ftp.ensembl.org/pub/release-${release}/gtf/mus_musculus/Mus_musculus.GRCm39.${release}.gtf.gz

# unzip files
gzip -d Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
gzip -d Mus_musculus.GRCm39.${release}.gtf.gz

# create a map of ENSEMBL gene IDs to gene names
awk 'FS="\t"{if($3=="transcript")print $9}' Mus_musculus.GRCm39.${release}.gtf | \
awk 'OFS="\t"{if($9=="gene_name") print $2, $6, $10; else print $2, $6, "NA"}' | \
sed 's/"//g; s/\;//g' > Mus_musculus.GRCm39.${release}.id_map.txt
