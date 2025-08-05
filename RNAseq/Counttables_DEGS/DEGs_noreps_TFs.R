# Install required packages (only run once)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomeInfoDbData")
BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))
install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))

# Load libraries
library(limma)
library(edgeR)
library(data.table)
library(pheatmap)
library(ggplot2)

# Set working directory
setwd("C:/Users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq")

#1. Read FeatureCounts Output
file_names <- list.files(path = ".", pattern = "featureCounts_gene\\.txt$", recursive = TRUE, full.names = TRUE)
# Read all count files into a DGEList
dgls <- readDGE(file_names, columns = c(1, 7), skip = 1)
counts_raw <- as.data.frame(dgls$counts)
# Rename sample columns
sample_names <- gsub("_featureCounts_gene", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sample_names

# 2. Clean and Filter Count Table
count_tbl <- counts_raw  # Simplified since we don't need to write/read from file
# Filter low-expressed genes (keep if expressed in ≥80% of samples)
perc_keep <- 0.8
gene_keep <- rowSums(count_tbl > 0) >= ceiling(perc_keep * ncol(count_tbl))
count_tbl_filtered <- count_tbl[gene_keep, ]

# 3, Create Sample Metadata
meta <- data.frame(SampleID = colnames(count_tbl),
                   CellLine = c("Control", "Y2", "YT"))
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

# 7. Create heatmap with specific genes of interest from CSV file
# Import genes of interest from a CSV file
genes_file <- read.csv("Top_100_TFs.csv")  # Replace with your actual file name
genes_of_interest <- genes_file$Gene_Symbol  # Adjust column name as needed

# Check which genes are present in your dataset
genes_present <- genes_of_interest[genes_of_interest %in% rownames(v$E)]
print(paste("Found", length(genes_present), "out of", length(genes_of_interest), "genes in the dataset"))

# If some genes are missing, print which ones
# if(length(genes_present) < length(genes_of_interest)) {
#   missing_genes <- genes_of_interest[!genes_of_interest %in% rownames(v$E)]
#   print("Missing genes:")
#   print(missing_genes)
# }

# Create a matrix with only those genes
log_counts_selected <- v$E[genes_present, ]
log_counts_selected_scaled <- t(scale(t(log_counts_selected)))

# Sample annotations for columns
annotation_col <- data.frame(CellLine = meta$CellLine)
rownames(annotation_col) <- meta$SampleID

# Generate the heatmap with just your genes of interest
mymap <- pheatmap(log_counts_selected_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         treeheight_row = 0,
         treeheight_col = 0,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         fontsize_row = 12,
         fontsize_col = 10,
         fontsize_main = 14,
         main = "Expression of Transcription Factors",
         filename = "Selected_genes_heatmap.png",
         width = 8,
         height = 6)

# Save using proper method for pheatmap objects
png("WT_Y2_YT2_100rnaseq_TFs.png", width = 10, height = 12, units = "in", res = 300)
print(mymap)
dev.off()
