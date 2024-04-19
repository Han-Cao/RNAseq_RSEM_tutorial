#!/bin/bash

#SBATCH -p q1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -o build_ref.log

###### SETTINGS ######
# set the input and output directories
ref_fa=/path/to/reference/genome
ref_gtf=/path/to/reference/annotation
outdir=/path/to/output/reference/
# set the prefix for the output reference files
prefix=rsem
# final reference to use in downstream analysis are:
# STAR: /path/to/output/reference/
# RSEM: /path/to/output/reference/rsem

# set number of threads to use
threads=40

# set reference
star_ref="/path/to/star_ref"
rsem_ref="/path/to/rsem_ref"

###### ANALYSIS ######
rsem-prepare-reference --gtf $ref_gtf \
                       --star \
                       -p $threads \
                       $ref_fa \
                       ${outdir}/${prefix}