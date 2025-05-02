# Load Required Libraries

library(DESeq2)
library(stringr)
library(ggplot2)
library(dplyr)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tibble)
library(AnnotationDbi)
library(ggrepel)

# Load raw count matrix

counts <- read.table("/Users/matthewng/Bioinformatics Spring 2025/HTP_Project_Files/GSE169687_raw_counts_GRCh38.p13_NCBI.tsv.gz",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)

# Load and Clean Metadata


# Read the series matrix as raw lines
lines <- readLines("/Users/matthewng/Bioinformatics Spring 2025/HTP_Project_Files/GSE169687_series_matrix(1).txt")

# Extract the annotation block (starts with sample titles)
start <- grep("!Sample_title", lines)
end   <- grep("!series_matrix_table_end", lines)
sample_lines <- lines[start:(end - 1)]

# Load the annotation block into a wide data frame
annotation_raw <- read.delim(text = sample_lines, header = FALSE, sep = "\t", stringsAsFactors = FALSE)

# Transpose and format the data
sample_data_t <- as.data.frame(t(annotation_raw[, -1]), stringsAsFactors = FALSE)
colnames(sample_data_t) <- annotation_raw[, 1]
sample_data_t$sample_id <- annotation_raw[annotation_raw$V1 == "!Sample_geo_accession", -1] %>% as.character()

# Identify the repeated characteristics columns
repeated_cols <- which(colnames(sample_data_t) == "!Sample_characteristics_ch1")

# Rename based on inspection order (confirmed via previous output)
colnames(sample_data_t)[repeated_cols] <- c(
  "age",
  "sex",
  "severity",        
  "rin",
  "sample_number",
  "tissue",
  "subject_id",
  "timepoint"        
)
sample_data_t <- as.data.frame(sample_data_t, stringsAsFactors = FALSE)


# Build cleaned metadata frame
meta_clean <- sample_data_t %>%
  dplyr::select(sample_id, timepoint, severity)%>%
  dplyr::mutate(
    sample_id = as.character(unlist(sample_id)),
    timepoint = str_replace(timepoint, "timepoint: ", ""),
    severity  = str_replace(severity, "disease severity: ", ""),
    condition = ifelse(severity == "Healthy", "Healthy", "Recovered"),
    timepoint = case_when(
      condition == "Healthy" ~ "Control",  
      TRUE ~ timepoint
    )
  ) %>%
  dplyr::filter(!is.na(sample_id), sample_id %in% colnames(counts), timepoint != "")


# Check Output

head(meta_clean)
table(meta_clean$timepoint)
table(meta_clean$condition)


################################################################################

# Filter metadata to include ALL Healthy and Recovered samples

meta_binary <- meta_clean %>%
  dplyr::filter(condition %in% c("Healthy", "Recovered")) %>%
  dplyr::mutate(
    sample_id = as.character(sample_id)
  )


# Match expression data to metadata

# Ensure the sample IDs are present in the count matrix
valid_samples <- intersect(meta_binary$sample_id, colnames(counts))
meta_binary <- meta_binary[meta_binary$sample_id %in% valid_samples, ]
counts_binary <- counts[, valid_samples]
rownames(meta_binary) <- meta_binary$sample_id


# Set factors for DESeq2

meta_binary$condition <- factor(meta_binary$condition, levels = c("Healthy", "Recovered"))


# Create DESeq2 object

dds_binary <- DESeqDataSetFromMatrix(
  countData = counts_binary,
  colData = meta_binary,
  design = ~ condition
)


# Run DESeq2 Analysis

dds_binary <- DESeq(dds_binary)


# Extract DEGs (Recovered vs Healthy)

res_binary <- results(dds_binary, contrast = c("condition", "Recovered", "Healthy"))
summary(res_binary)


res_binary_ordered <- res_binary[order(res_binary$pvalue), ]
write.csv(as.data.frame(res_binary_ordered), "DEGs_Recovered_vs_Healthy_ALL.csv", row.names = TRUE)

# Visualization (MA Plot)
plotMA(res_binary, ylim = c(-2, 2), main = "Recovered vs Healthy (All Timepoints)")

# Preview Top Genes

head(res_binary_ordered, 10)



################################################################################

# Convert DESeq2 result to data frame with gene IDs
res_df <- as.data.frame(res_binary) %>%
  rownames_to_column("ensembl")

# Remove version suffix from Ensembl IDs 
res_df$ensembl_clean <- gsub("\\..*", "", res_df$ensembl)

# Map to gene symbols
res_df$symbol <- mapIds(
  org.Hs.eg.db,
  keys = res_df$ensembl_clean,
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Rearrange and filter if needed
res_df <- res_df %>%
  dplyr::select(ensembl, symbol, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)

# Save with symbols
write.csv(res_df, "DEGs_with_gene_symbols.csv", row.names = FALSE)

# View significant ones with gene symbols
res_sig <- res_df %>% filter(padj < 0.05, abs(log2FoldChange) > 1)
write.csv(res_sig, "Significant_DEGs_with_gene_symbols.csv", row.names = FALSE)
head(res_sig)

vsd <- vst(dds_binary, blind = FALSE)
vsd_mat <- assay(vsd)  

# Check a few gene IDs
head(rownames(vsd_mat))

# Map Entrez IDs to gene symbols (may result in NA)
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(vsd_mat),
  column = "SYMBOL",
  keytype = "ENTREZID",
  multiVals = "first"
)

# Replace NAs with original Entrez ID
gene_symbols_filled <- ifelse(is.na(gene_symbols), rownames(vsd_mat), gene_symbols)

# Assign to rownames
rownames(vsd_mat) <- gene_symbols_filled



# Subset only to significant DEGs

degs <- read.csv("Significant_DEGs_with_gene_symbols.csv")
sig_genes <- degs$symbol

vsd_deg_only <- vsd_mat[rownames(vsd_mat) %in% sig_genes, ]


# Transpose for ML and add labels

vsd_ml <- t(vsd_deg_only)

# Label: 0 = Healthy, 1 = Recovered
labels <- ifelse(meta_binary$condition == "Recovered", 1, 0)

# Create dataframe and attach labels
ml_df <- as.data.frame(vsd_ml)
ml_df$label <- labels[match(rownames(ml_df), meta_binary$sample_id)]

# Save for ML
write.csv(ml_df, "ML_input_expression_symbols.csv", row.names = TRUE)

################################################################################

# Load DEG results
degs <- read.csv("DEGs_Recovered_vs_Healthy_ALL.csv", row.names = 1)

# Add -log10 adjusted p-value column
degs$log10_padj <- -log10(degs$padj)

# Define significance thresholds
logFC_cutoff <- 0.5
padj_cutoff <- 0.05

# Annotate genes by significance
degs$category <- "Not significant"
degs$category[degs$log2FoldChange >= logFC_cutoff & degs$padj < padj_cutoff] <- "Upregulated"
degs$category[degs$log2FoldChange <= -logFC_cutoff & degs$padj < padj_cutoff] <- "Downregulated"

# Keep top labels only (top N by significance)
degs$label <- ifelse(
  degs$padj < padj_cutoff & abs(degs$log2FoldChange) >= logFC_cutoff,
  rownames(degs),
  NA
)


degs <- read.csv("DEGs_Recovered_vs_Healthy_ALL.csv", row.names = 1)
degs <- tibble::rownames_to_column(degs, var = "ensembl")

res_df <- read.csv("DEGs_with_gene_symbols.csv")  # this file contains ensembl + symbol
res_df$ensembl <- as.character(res_df$ensembl)

degs_merged <- left_join(degs, res_df[, c("ensembl", "symbol")], by = "ensembl")

logFC_cutoff <- 0.5
padj_cutoff <- 0.05

degs_merged$log10_padj <- -log10(degs_merged$padj)

degs_merged$category <- "Not significant"
degs_merged$category[degs_merged$log2FoldChange >= logFC_cutoff & degs_merged$padj < padj_cutoff] <- "Upregulated"
degs_merged$category[degs_merged$log2FoldChange <= -logFC_cutoff & degs_merged$padj < padj_cutoff] <- "Downregulated"

degs_merged$label <- ifelse(
  degs_merged$padj < padj_cutoff & abs(degs_merged$log2FoldChange) >= logFC_cutoff,
  degs_merged$symbol,
  NA
)


# Plot
ggplot(degs_merged, aes(x = log2FoldChange, y = log10_padj, color = category)) +
  geom_point(alpha = 0.8, size = 1.8) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not significant" = "gray")) +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20) +
  labs(
    title = "Volcano Plot: Recovered vs Healthy",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    color = "Significance"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "top"
  )
table(degs_merged$category)


colnames(degs)
colnames(ml_df)
################################################################################

# Use all DEGs from volcano plot
degs_all <- degs_merged %>%
  filter(category %in% c("Upregulated", "Downregulated"))

write.csv(degs_all$symbol, "all_DEGs_symbols.csv", row.names = FALSE)

# Don't filter out NA gene symbols
gene_symbols_filled <- ifelse(is.na(gene_symbols), rownames(vsd_mat), gene_symbols)
rownames(vsd_mat) <- gene_symbols_filled

# Keep all genes including those with unmapped symbols
vsd_deg_all <- vsd_mat[rownames(vsd_mat) %in% degs_all$symbol, ]

vsd_ml_all <- t(vsd_deg_all)
ml_df_all <- as.data.frame(vsd_ml_all)
ml_df_all$label <- labels[match(rownames(ml_df_all), meta_binary$sample_id)]

# Use this for ML models
write.csv(ml_df_all, "ML_input_all_volcano_genes.csv", row.names = TRUE)


