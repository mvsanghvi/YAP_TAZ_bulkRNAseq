#Set Up environment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))

library(edgeR)
library(limma)
library(data.table)
library(pheatmap)
library(ggplot2)

setwd("C:/Users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq")

#1. Read FeatureCounts Output
file_names <- list.files(path = ".", pattern = "featureCounts_gene\\.txt$", recursive = TRUE, full.names = TRUE)

# Read all count files into a DGEList
dgls <- readDGE(file_names, columns = c(1, 7), skip = 1)
counts_raw <- as.data.frame(dgls$counts)

# Rename sample columns
sample_names <- gsub("_featureCounts_gene", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sample_names

# Save raw counts (optional)
fwrite(counts_raw, "counts_raw.tsv", sep = "\t", row.names = TRUE)

# 2. Clean and Filter Count Table

count_tbl <- fread("counts_raw.tsv", data.table = FALSE)
rownames(count_tbl) <- count_tbl[[1]]
count_tbl <- count_tbl[, -1]

# Filter low-expressed genes (keep if expressed in ≥80% of samples)
perc_keep <- 0.8
gene_keep <- rowSums(count_tbl > 0) >= ceiling(perc_keep * ncol(count_tbl))
count_tbl_filtered <- count_tbl[gene_keep, ]

# 3, Create Sample Metadata
meta <- data.frame(SampleID = colnames(count_tbl),
                   CellLine = c("Control", "Y2", "YT"))  # Edit as needed
rownames(meta) <- meta$SampleID

# 4. Create DGEList and Normalize
dge <- DGEList(counts = count_tbl_filtered, samples = meta)
dge <- calcNormFactors(dge)

# 5. Create Design Matrix for voom
design <- model.matrix(~ 0 + dge$samples$CellLine)
colnames(design) <- levels(factor(dge$samples$CellLine))

# 6. voom Transformation (No DE)
v <- voom(dge, design, plot = TRUE)
v$targets <- meta  # Attach metadata for easy access

# Identify genes with Ensembl IDs (ENSG pattern)
ensembl_pattern <- grep("^ENSG", rownames(v$E), value = TRUE)

# Create a filtered expression matrix without Ensembl-only genes
v_filtered <- v
v_filtered$E <- v$E[!rownames(v$E) %in% ensembl_pattern, ]

# Now calculate variance using the filtered data
gene_vars <- apply(v_filtered$E, 1, var)

# Select top 100 most variable genes (now only from genes with proper symbols)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]
log_counts_top <- v_filtered$E[top_genes, ]

# Scale for heatmap (row-wise)
log_counts_top_scaled <- t(scale(t(log_counts_top)))

# Sample annotations for columns
annotation_col <- data.frame(CellLine = meta$CellLine)
rownames(annotation_col) <- meta$SampleID

# Create the heatmap and save it to a variable
mymap <- pheatmap(log_counts_top_scaled,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  treeheight_row = 0,
                  treeheight_col = 0,
                  annotation_col = annotation_col,
                  show_rownames = TRUE,
                  fontsize_row = 7,
                  fontsize_col = 10,
                  fontsize_main = 14,
                  main = "Top 100 Most Variable Genes (Named Genes Only)")

# Save using proper method for pheatmap objects
png("WT_Y2_YT2_100rnaseq_noENS.png", width = 10, height = 12, units = "in", res = 300)
print(mymap)
dev.off()
