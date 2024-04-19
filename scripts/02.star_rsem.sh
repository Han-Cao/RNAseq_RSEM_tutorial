#!/bin/bash
#SBATCH -p q1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 40
#SBATCH -o star_rsem.log

# star 2.7.10b
# rsem 1.3.3

[[ ! -d "${outpath}/${sample}" ]] && mkdir -p "${outpath}/${sample}"

###### SETTINGS ######
# set the input and output directories
inpath=/path/to/fastp/output/
outpath=/path/to/star_rsem/output/
# set number of threads to use
threads=40

# set reference
star_ref="/path/to/star_ref"
rsem_ref="/path/to/rsem_ref"

###### ANALYSIS ######
# loop all samples in the input directory
ls ${inpath} | while read sample;
do
    # create the output directory if it does not exist
    [[ ! -d "${outpath}/${sample}" ]] && mkdir -p "${outpath}/${sample}"

    # extract all input fastq names
    read1=$(ls ${inpath}/${sample}/*_1_paired.fq.gz | paste -s -d ',')
    read2=$(echo $read1 | sed 's/_1_paired.fq.gz/_2_paired.fq.gz/g')

    # STAR
    STAR --genomeDir $star_ref \
    --outSAMunmapped Within  --outFilterType BySJout  --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20  --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverLmax 0.04  --alignIntronMin 20  --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000  --alignSJoverhangMin 8  --alignSJDBoverhangMin 1 \
    --sjdbScore 1  --runThreadN $threads  --genomeLoad NoSharedMemory  --outSAMtype BAM Unsorted \
    --quantMode TranscriptomeSAM  --outSAMheaderHD @HD VN:1.4 SO:unsorted \
    --outFileNamePrefix "${outpath}/${sample}/${sample}.STAR."  --readFilesCommand zcat \
    --readFilesIn $read1 $read2

    # sort genome bam
    samtools sort -@ $threads -O BAM --write-index \
    -o "${outpath}/${sample}/${sample}.STAR.Aligned.toGenome.sort.bam" \
    "${outpath}/${sample}/${sample}.STAR.Aligned.out.bam"

    # remove unsorted genome bam file
    if [[ -e "${outpath}/${sample}/${sample}.STAR.Aligned.out.bam" ]]; then
        rm "${outpath}/${sample}/${sample}.STAR.Aligned.out.bam"
    fi

    # rsem
    rsem-calculate-expression --alignments --no-bam-output -p $threads \
    --paired-end --forward-prob 0.5 \
    "${outpath}/${sample}/${sample}.STAR.Aligned.toTranscriptome.out.bam" \
    $rsem_ref \
    "${outpath}/${sample}/${sample}.rsem"
done

