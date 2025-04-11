#Set Up environment 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))

library(limma)
library(edgeR)
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
v$targets <- meta

# Install and use ReactomePA for pathway analysis
if(interactive()) {
  cat("\nWould you like to try installing ReactomePA for pathway analysis? (y/n): ")
  answer <- readline()
  if(tolower(answer) == "y") {
    if(!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
    }
    BiocManager::install("ReactomePA")
    cat("Please run the ReactomePA code separately after installation completes.\n")
  }
}

# Load necessary libraries
library(ReactomePA)
library(biomaRt)
library(enrichplot)  # For visualization
library(ggplot2)
library(pheatmap)

# 7. Pathway Analysis
# Create a vector of gene variances
gene_vars <- apply(v$E, 1, var)
names(gene_vars) <- rownames(v$E)

# Sort genes by variance
ranked_genes <- sort(gene_vars, decreasing = TRUE)
top_genes <- names(ranked_genes)[1:500]  # Top 500 most variable genes

# Convert gene symbols to Entrez IDs using biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")  # For human (change for other organisms)
entrez_ids <- getBM(attributes = c("external_gene_name", "entrezgene_id"),
                    filters = "external_gene_name",
                    values = top_genes,
                    mart = mart)

# Remove any NA values
entrez_ids <- entrez_ids[!is.na(entrez_ids$entrezgene_id),]

# Run Reactome pathway analysis
reactome_results <- enrichPathway(gene = entrez_ids$entrezgene_id,
                                  pvalueCutoff = 0.05,
                                  readable = TRUE)

# Print top pathways
print(head(reactome_results, 10))

# Visualize Reactome pathway results
# Dotplot showing enriched pathways
pdf("Reactome_enrichment_dotplot.pdf", width = 10, height = 8)
dotplot(reactome_results, showCategory = 15, 
        title = "Reactome Pathway Enrichment")
dev.off()

# Create pathway-gene heatmap
# Extract top pathways and their genes
top_pathways <- head(reactome_results@result$Description, 10)  # Top 10 pathways
pathway_genes <- list()

# Get genes for each pathway
for(pathway in top_pathways) {
  pathway_idx <- which(reactome_results@result$Description == pathway)
  genes_in_pathway <- strsplit(reactome_results@result$geneID[pathway_idx], "/")[[1]]
  # Map back to gene symbols
  gene_symbols <- entrez_ids$external_gene_name[match(genes_in_pathway, entrez_ids$entrezgene_id)]
  pathway_genes[[pathway]] <- gene_symbols
}

# Create a matrix for heatmap (pathways x genes)
# Get all unique genes in the top pathways
all_pathway_genes <- unique(unlist(pathway_genes))
all_pathway_genes <- all_