#Set Up environment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db"))

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

# 6. voom Transformation (No DEGs)
v <- voom(dge, design, plot = TRUE)
v$targets <- meta  # Attach metadata for easy access
fwrite(v$targets, "Voom_YT.tsv", sep = "\t", row.names = T)

# Export the voom-transformed expression values
voom_expr <- as.data.frame(v$E)
voom_expr$gene_id <- rownames(voom_expr)  # Add gene IDs as a column
fwrite(voom_expr, "voom_expression_data.tsv", sep="\t", row.names=FALSE)

# 7. Heatmap of Top Variable Genes
# Calculate variance of each gene
gene_vars <- apply(v$E, 1, var)
gene_var_df <- data.frame(gene_id = names(gene_vars), variance = gene_vars)
gene_var_df <- gene_var_df[order(-gene_var_df$variance),]  # Sort by variance
fwrite(gene_var_df, "gene_variance.tsv", sep="\t", row.names=FALSE)

# Select top 100 most variable genes
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]
log_counts_top <- v$E[top_genes, ]

# Scale for heatmap (row-wise)
log_counts_top_scaled <- t(scale(t(log_counts_top)))

# Sample annotations for columns
annotation_col <- data.frame(CellLine = meta$CellLine)
rownames(annotation_col) <- meta$SampleID

# Plot heatmap
mymap <- pheatmap(log_counts_top_scaled,
                  cluster_rows = TRUE,
                  cluster_cols = TRUE,
                  treeheight_row = 0,  # Set dendrogram height to 0 (hides it)
                  treeheight_col = 0,  # Set dendrogram height to 0 (hides it)
                  annotation_col = annotation_col,
                  show_rownames = TRUE,
                  fontsize_row = 5,
                  fontsize_col = 5,
                  main = "Top 100 Most Variable Genes")

# Save to image
ggsave("WT_Y2_YT2_rnaseq.png", plot=mymap)

# save_heatmap_pdf <- function(filename, mat, annotation) {
#   pdf(filename, width = 8, height = 10)
#   pheatmap(mat,
#            cluster_rows = TRUE,
#            cluster_cols = TRUE,
#            annotation_col = annotation,
#            show_rownames = FALSE,
#            fontsize_col = 10,
#            main = "Top 100 Most Variable Genes")
#   dev.off()
# }
# 
# save_heatmap_pdf("heatmap_top100.pdf", log_counts_top_scaled, annotation_col)

#Export "100 most Variable Genes" Genes to CSV
# After you've identified the top 100 most variable genes
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]

# Create a data frame with gene IDs and their variance
top_genes_df <- data.frame(
  Gene_ID = top_genes,
  Variance = gene_vars[top_genes]
)

# If you've already mapped Ensembl IDs to gene symbols, add those too
if(exists("gene_map")) {
  top_genes_df$Gene_Symbol <- sapply(top_genes_df$Gene_ID, function(id) {
    if(id %in% names(gene_map) && gene_map[id] != "") {
      return(gene_map[id])
    } else {
      return(id)  # Keep original ID if no symbol available
    }
  })
}

# Add expression values for each sample
top_genes_df <- cbind(top_genes_df, v$E[top_genes, ])

# Export as CSV
write.csv(top_genes_df, "top_100_variable_genes.csv", row.names = FALSE)
