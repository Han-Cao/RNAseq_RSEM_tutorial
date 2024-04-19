#!/bin/bash
#SBATCH -p q1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -o fastp.log

# fastp 0.23.4

###### SETTINGS ######
# set the input and output directories
inpath=/path/to/fastq/input/
outpath=/path/to/fastp/output/
# set number of threads to use
# the maximum number of threads for fastp is 16
threads=16

###### ANALYSIS ######
# loop all samples in the input directory
ls ${inpath} | while read sample;
do
    # create the output directory if it does not exist
    [[ ! -d "${outpath}/${sample}" ]] && mkdir -p "${outpath}/${sample}"
    # loop all fastq files in the sample directory
    ls ${inpath}${sample}/*_1.fq.gz | while read fastq
    do
        # set names for fq1, fq2, and output name
        fastq=${fastq%_1.fq.gz}
        fq1=${fastq}_1.fq.gz
        fq2=${fastq}_2.fq.gz
        name=$(basename $fastq)

        # run fastp to preprocess the fastq files
        fastp --thread $threads \
        --in1 $fq1 --in2 $fq2 \
        --out1 "${outpath}/${sample}/${name}_1_paired.fq.gz" \
        --out2 "${outpath}/${sample}/${name}_2_paired.fq.gz" \
        --unpaired1 "${outpath}/${sample}/${name}_1_unpaired.fq.gz" \
        --unpaired2 "${outpath}/${sample}/${name}_2_unpaired.fq.gz" \
        --length_required 40 \
        --verbose \
        --json "${outpath}/${sample}/${name}_fastp.json" \
        --html "${outpath}/${sample}/${name}_fastp.html"
    done
done
