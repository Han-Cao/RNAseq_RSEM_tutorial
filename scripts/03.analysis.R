library(readr)
library(dplyr)
library(tximport)
library(apeglm)
library(DESeq2)

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
write_tsv(df_res, "output/DEG_result_lfcShrink.txt")

