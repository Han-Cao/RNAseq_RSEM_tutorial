# RNAseq data analysis

## Table of contents
- [Introduction](#introduction)
- [Quick start](#quick-start)
- [Prepare for analysis](#prepare-for-analysis)
  - [Setup environment](#setup-environment)
  - [Work with SLURM](#work-with-slurm)
  - [Raw data requirements](#raw-data-requirements)
  - [Build reference](#build-reference)
- [Data analysis](#data-analysis)
    - [Data preprocessing](#data-preprocessing)
    - [Align RNAseq data and quantify gene expression](#align-rnaseq-data-and-quantify-gene-expression)
    - [Differential gene expression analysis](#differential-gene-expression-analysis)


## Introduction
This repository contains the scripts for RNAseq data processing and differential gene expression analysis. You can find all the scripts to use in the `scripts` folder. The `example` folder contains the example metadata that is required to run differential gene expression analysis.

## Quick start
When first run the analysis, please do:
- [Setup environment](#setup-environment)
- [Build reference](#build-reference)

After activate the conda enviroment, you can sequentially run the following scripts to quantify gene expression and perform differential gene expression analysis:
- [Data preprocessing](#data-preprocessing) : [01.fastp.sh](scripts/01.fastp.sh)
- [Align RNAseq data and quantify gene expression](#align-rnaseq-data-and-quantify-gene-expression): [02.star_rsem.sh](scripts/02.star_rsem.sh)
- [Differential gene expression analysis](#differential-gene-expression-analysis): [03.analysis.R](scripts/03.analysis.R)

If working on a cluster with SLURM, please always submit jobs with `sbatch your_script.sh` command.

## Prepare for analysis
### Setup environment
Before you start, please make sure conda is installed on the system. If not, please follow the instructions  to install [conda](https://docs.anaconda.com/free/miniconda/miniconda-install/) or [mamba](https://mamba.readthedocs.io/en/latest/installation/micromamba-installation.html).

Once conda has been installed and initialized, you can install the required tools in a new conda environment with the following command:
```bash
conda create -n rnaseq -c bioconda -c conda-forge fastp=0.23.4 star=2.7.10b rsem=1.3.3 samtools=1.20
```
Please remember to activate the environment before any analysis:
```bash
conda activate rnaseq
```

### Work with SLURM
If working on a cluster with SLURM, please **always** submit jobs with the following command:
```bash
sbatch script.sh
```
You need to specify the resources to use in the `script.sh` file. For example:
```bash
#!/bin/bash
#SBATCH -p q1
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -o output.log
```
This means submitting the job to the `q1` partition, requesting 1 node for 1 task using 20 cores. The log will be saved in the `output.log` file. Please do not request more resource than the script need and the capability of the server.

### Raw data requirements
The scripts are designed to analyze **paired-end** RNAseq data with the following data structure:
- All fastq files are ended with `_1.fq.gz` and `_2.fq.gz`
- All fastq files from the same sample are stored in the same folder.
- The parent folder of all samples **ONLY** contains the sample folders


Correct:

```bash
- /path/to/fastq/input/
    - Sample1
        - Sample1_1.fq.gz
        - Sample1_2.fq.gz
    - Sample2
        - Sample2_1.fq.gz
        - Sample2_2.fq.gz
```

Incorrect:

```bash
- /path/to/fastq/input/
    - Sample1
        - Sample1_R1.fq.gz  # _R1.fq.gz not works
        - Sample1_R2.fq.gz  # _R2.fq.gz not works
    - Sample2
        - Sample2_1.fq.gz
        - Sample2_2.fq.gz
    - Sample3_1.fq.gz       # Create a Sample3 folder for this
    - Sample3_2.fq.gz       # Create a Sample3 folder for this
    - anyfiles.txt          # any other files will cause error
```

If your data is in different structure, please rename the files or modify the scripts accordingly.


### Build reference
Use the [download_files.sh](scripts/download_files.sh) to download the reference genome and annotation files for Mus musculus from the Ensembl FTP website. The script will also build a map of Ensembl gene IDs to gene names. You can modify the script to download files for other species. 

```bash
#!/bin/bash

###### SETTINGS ######
# ENSEMBL release version
# change if you want to use a different version
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
```

The above script will generate following reference files:
- reference genome: `Mus_musculus.GRCm39.dna.primary_assembly.fa`
- transcript annotation: `Mus_musculus.GRCm39.111.gtf`
- ID map between ENSEMBL ID and gene symbol: `Mus_musculus.GRCm39.111.id_map.txt`

After downloading the files, the STAR index and RSEM reference can be built by using the [00.build_ref.sh](scripts/00.build_ref.sh) script. Please remember to modify the input and output directories in the script before running. Please use the following path as reference input for STAR and RSEM
- STAR: `/path/to/output/reference/`
- RSEM: `/path/to/output/reference/rsem`

```bash
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
```

## Data analysis
### Data preprocessing
Run the [01.fastp.sh](scripts/01.fastp.sh) script to preprocess the raw fastq files and perform quality control. The script will automatically trim adapters and filter low-quality reads.

Clean reads will be saved in the following fastq files and should be used in the downstream analysis:
- `${outpath}/${sample}/${name}_1_paired.fq.gz`
- `${outpath}/${sample}/${name}_2_paired.fq.gz`

Quality control reports will be saved in:
- `${outpath}/${sample}/${name}_fastp.html` for browser view. An example can be found [here](https://htmlpreview.github.io/?https://github.com/Han-Cao/RNAseq_RSEM_tutorial/blob/master/example/fastp.html).
- `${outpath}/${sample}/${name}_fastp.json` for machine-readable format

```bash
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
```

### Align RNAseq data and quantify gene expression
Run the [02.star_rsem.sh](scripts/02.star_rsem.sh) script to align the preprocessed fastq files (`*paired.fq.gz`) to the reference genome and quantify gene expression. The script will generate the expression matrix for differential gene expression analysis:

- `${outpath}/${sample}/${sample}.rsem.genes.results` for gene-level expression
- `${outpath}/${sample}/${sample}.rsem.isoforms.results` for isoform-level expression

```bash
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
```

### Differential gene expression analysis
The [03.analysis.R](scripts/03.analysis.R) script performs differential gene expression analysis using the R package `DESeq2`. The following R packages are required:

```R
library(readr)
library(dplyr)
library(tximport)
library(apeglm)
library(DESeq2)
```

Before the analysis, you need to prepare a metadata file of study design, as well as a gene ID map to match ENSEMBL ID to gene symbol. The example files can be found in the `example` folder:

metadata file:
- The column names must be "sample", "file", and "condition" if you don't want to modify the R script
- If there are multiple variables to test, can add more columns and modify the design formula in the R script (see below)
- The values in the "condition" column need to match the values in `contrast` (see below)

| sample | file                          | condition  |
|--------|-------------------------------|------------|
| A      | /path/to/A.rsem.genes.results | Control    |
| B      | /path/to/B.rsem.genes.results | Control    |
| C      | /path/to/C.rsem.genes.results | Control    |
| D      | /path/to/D.rsem.genes.results | Treat      |
| E      | /path/to/E.rsem.genes.results | Treat      |
| F      | /path/to/F.rsem.genes.results | Treat      |

Gene ID map:
- This is generated in the step (Build reference)[#build-reference]
- Columns are ENSEMBL gene id, transcript id, and gene symbol **without** header

|                 |                 |          |
|-----------------|-----------------|----------|
| ENSG00000160072 | ENST00000673477 | ATAD3B   |
| ENSG00000160072 | ENST00000472194 | ATAD3B   |
| ENSG00000160072 | ENST00000378736 | ATAD3B   |
| ENSG00000160072 | ENST00000485748 | ATAD3B   |
| ENSG00000160072 | ENST00000474481 | ATAD3B   |
| ENSG00000160072 | ENST00000308647 | ATAD3B   |
| ENSG00000279928 | ENST00000624431 | DDX11L17 |

Then, then data can be loaded into R:
```R
###### SETTINGS ######
# set the input and output directories
sample_tab <- "sample_table.txt"
id_map <- "Mus_musculus.GRCm39.111.id_map.txt"

###### ANALYSIS ######
# read tables for samples and genes
df_sample <- read_tsv(sample_tab)
df_map_isoform <- read_tsv(id_map, col_names = FALSE)
colnames(df_map_isoform) <- c("Gene", "Isoform", "Symbol")
df_map_gene <- dplyr::select(df_map_isoform, Gene, Symbol) %>% dplyr::filter(!duplicated(Gene)) # remove duplicated records at gene level

# load gene expression matrix
rsem_files <- df_sample$file
names(rsem_files) <- df_sample$sample
txi <- tximport(rsem_files, type = "rsem", txIn = FALSE, txOut = FALSE)
# technically, RSEM produce effect gene length of 0 when gene length < read length
# it can be converted to 1 to allow DESeq2 estimate library size
txi$length[txi$length == 0] <- 1
```

The gene expression matrix can be extracted into a table:
- gene expression matrix in counts: `output/gene_expression_matrix_counts.txt`
- gene expression matrix in TPM: `output/gene_expression_matrix_tpm.txt`

```R
# convert expression matrix to table and map gene ID
convert_expr_mat <- function(x, type, df_map){
  if(type == "Gene"){
    df_new <- as_tibble(x, rownames="Gene")
    df_new <- right_join(df_map, df_new, by = "Gene")
  }
  else if (type=="Isoform") {
    df_new <- as_tibble(x, rownames="Isoform")
    df_new <- right_join(df_map, df_new, by="Isoform")
  }
  return(df_new)
}
df_rsem_gene_counts <- convert_expr_mat(txi$counts, "Gene", df_map_gene)
df_rsem_gene_tpm <- convert_expr_mat(txi$abundance, "Gene", df_map_gene)

write_tsv(df_rsem_gene_counts, "output/gene_expression_matrix_counts.txt")
write_tsv(df_rsem_gene_tpm, "output/gene_expression_matrix_tpm.txt")
```

The differential gene expression analysis can be performed by below commands. You may need to modify some parameters to match the study design:
- `contract` should be `c("column_to_test", "test_condition", "control_condition")`
- `coef` should be `"[column_to_test]_[test_condition]_vs_[control_condition]"`

```R
# DEG analysis
# load txi to DESeq2
dds <- DESeqDataSetFromTximport(txi, colData = df_sample, design = ~ condition)
# run analysis
dds <- DESeq(dds)
# DEG results
# the contrast is defined as column_name_in_table, test_condidtion, control_condition
# these values MUST be identical to the values in the sample_table.txt
contrast <- c("condition", "Treat", "Control")
res <- results(dds, contrast = contrast, alpha = 0.05)
resLFC <- lfcShrink(dds, coef = "condition_Treat_vs_Control", type = "apeglm")

# convert the results to table and map gene names
df_res <- as_tibble(res, rownames = "Gene")
df_resLFC <- as_tibble(resLFC, rownames = "Gene")

df_res <- right_join(df_map_gene, df_res, by=c("Gene"))
df_resLFC <- right_join(df_map_gene, df_resLFC, by=c("Gene"))

write_tsv(df_res, "output/DEG_result.txt")
write_tsv(df_resLFC, "output/DEG_result_lfcShrink.txt")
```

Finally, the differential gene expression results are saved in:
- `output/DEG_result.txt` for original DESeq2 results
- `output/DEG_result_lfcShrink.txt` for DESeq2 results with log fold change shrinkage. This is useful for visualization and ranking genes, especially for lowly-expressed genes with high variability.
