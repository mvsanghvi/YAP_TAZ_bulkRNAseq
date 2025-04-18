#Set Up environment 
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db")
BiocManager::install("GenomeInfoDbData")
BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db",
                       "ReactomePA", "DOSE", "pathview", "enrichplot"))
install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse", "ComplexHeatmap"))
library(limma)
library(edgeR)
library(data.table)
library(pheatmap)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(ComplexHeatmap)
library(RColorBrewer)
library(biomaRt)
library(tidyverse)

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

# 3. Create Sample Metadata
meta <- data.frame(SampleID = colnames(count_tbl),
                   CellLine = c("Control", "Y2", "YT"))  # Edit as needed
rownames(meta) <- meta$SampleID

# 4. Create DGEList and Normalize
dge <- DGEList(counts = count_tbl_filtered, samples = meta)
dge <- calcNormFactors(dge)

# 5. Create Design Matrix for voom
design <- model.matrix(~ 0 + dge$samples$CellLine)
colnames(design) <- levels(factor(dge$samples$CellLine))

# 6. voom Transformation
v <- voom(dge, design, plot = TRUE)
v$targets <- meta

# Instead of contrast-based DE analysis, use a simpler approach for unreplicated data

# Get normalized expression data
norm_expr <- v$E

# Get normalized expression data
norm_expr <- v$E

# Make sure column names are correctly aligned with metadata
print(colnames(norm_expr))  # Check what your actual column names are
print(v$targets$SampleID)   # Check your sample IDs
print(v$targets$CellLine)   # Check your cell line designations

# Extract expression based on actual column names and cell line info
control_samples <- v$targets$SampleID[v$targets$CellLine == "Control"]
Y2_samples <- v$targets$SampleID[v$targets$CellLine == "Y2"]
YT_samples <- v$targets$SampleID[v$targets$CellLine == "YT"]

# Calculate fold changes manually between conditions
# For each cell line, we need to use the correct sample name(s)
log2FC_Y2vsControl <- norm_expr[, Y2_samples] - norm_expr[, control_samples]
log2FC_YTvsControl <- norm_expr[, YT_samples] - norm_expr[, control_samples]
log2FC_Y2vsYT <- norm_expr[, Y2_samples] - norm_expr[, YT_samples]

# Create result tables with fold changes (but no p-values)
Y2vsControl_DE <- data.frame(
  logFC = log2FC_Y2vsControl,
  gene_id = rownames(norm_expr),
  row.names = rownames(norm_expr)
)

YTvsControl_DE <- data.frame(
  logFC = log2FC_YTvsControl,
  gene_id = rownames(norm_expr),
  row.names = rownames(norm_expr)
)

Y2vsYT_DE <- data.frame(
  logFC = log2FC_Y2vsYT,
  gene_id = rownames(norm_expr),
  row.names = rownames(norm_expr)
)

# Sort by absolute fold change
Y2vsControl_DE <- Y2vsControl_DE[order(abs(Y2vsControl_DE$logFC), decreasing = TRUE),]
YTvsControl_DE <- YTvsControl_DE[order(abs(YTvsControl_DE$logFC), decreasing = TRUE),]
Y2vsYT_DE <- Y2vsYT_DE[order(abs(Y2vsYT_DE$logFC), decreasing = TRUE),]

# Add a column for fold change direction
Y2vsControl_DE$direction <- ifelse(Y2vsControl_DE$logFC > 0, "up", "down")
YTvsControl_DE$direction <- ifelse(YTvsControl_DE$logFC > 0, "up", "down")
Y2vsYT_DE$direction <- ifelse(Y2vsYT_DE$logFC > 0, "up", "down")

# Save results
write.csv(Y2vsControl_DE, "Y2vsControl_DE.csv")
write.csv(YTvsControl_DE, "YTvsControl_DE.csv")
write.csv(Y2vsYT_DE, "Y2vsYT_DE.csv")

# Prepare gene lists based on fold change magnitude
prepare_gene_list <- function(de_results, id_map, lfc_threshold = 1) {

  # Sort genes by fold change
  de_sorted <- de_results[order(de_results$logFC, decreasing = TRUE),]
  
  # Extract top genes with large fold changes
  up_genes <- rownames(de_sorted[de_sorted$logFC >= lfc_threshold,])
  down_genes <- rownames(de_sorted[de_sorted$logFC <= -lfc_threshold,])
  
  # Map to Entrez IDs
  up_entrez <- id_map$entrezgene_id[id_map$ensembl_gene_id %in% up_genes]
  down_entrez <- id_map$entrezgene_id[id_map$ensembl_gene_id %in% down_genes]
  
  # Create ranked gene list for GSEA
  all_genes <- de_results$logFC
  names(all_genes) <- rownames(de_results)
  
  # Map to Entrez IDs, keeping logFC values
  entrez_map <- id_map[match(names(all_genes), id_map$ensembl_gene_id),]
  ranked_list <- all_genes[!is.na(entrez_map$entrezgene_id)]
  names(ranked_list) <- entrez_map$entrezgene_id[!is.na(entrez_map$entrezgene_id)]
  ranked_list <- sort(ranked_list, decreasing = TRUE)
  
  return(list(
    up_genes = up_entrez, 
    down_genes = down_entrez, 
    ranked_list = ranked_list
  ))
}

# Modify GO analysis function to analyze up and down regulated genes separately
run_go_analysis <- function(gene_list, ont = "BP", gene_universe = NULL) {
  # Up-regulated genes
  ego_up <- enrichGO(
    gene = gene_list$up_genes,
    universe = gene_universe,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  # Down-regulated genes
  ego_down <- enrichGO(
    gene = gene_list$down_genes,
    universe = gene_universe,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.1,
    readable = TRUE
  )
  
  return(list(up = ego_up, down = ego_down))
}

# Create heatmaps after fold change calculation

# Libraries needed for heatmaps (make sure these are loaded)
library(pheatmap)
library(RColorBrewer)

# 1. Heatmap of top differentially expressed genes by fold change
# Get top genes based on fold change magnitude
get_top_fc_genes <- function(fc_values, n = 50) {
  names(fc_values)[order(abs(fc_values), decreasing = TRUE)][1:min(n, length(fc_values))]
}

# Get top genes for each comparison
top_Y2vsControl <- get_top_fc_genes(log2FC_Y2vsControl)
top_YTvsControl <- get_top_fc_genes(log2FC_YTvsControl)
top_Y2vsYT <- get_top_fc_genes(log2FC_Y2vsYT)

# Combine top genes from all comparisons (unique)
all_top_genes <- unique(c(top_Y2vsControl, top_YTvsControl, top_Y2vsYT))

# Get expression data for top DE genes
de_expr_matrix <- norm_expr[all_top_genes, ]

# Try to add gene symbols if you have them (from biomaRt)
# If you don't have gene symbols yet, run this code:
library(biomaRt)
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_info <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = all_top_genes,
  mart = ensembl
)

# Create a mapping and rename rows
gene_symbols <- gene_info$external_gene_name
names(gene_symbols) <- gene_info$ensembl_gene_id
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

# Replace gene IDs with symbols where available
row_labels <- rownames(de_expr_matrix)
for (i in 1:length(row_labels)) {
  gene_id <- row_labels[i]
  if (gene_id %in% names(gene_symbols)) {
    row_labels[i] <- paste0(gene_symbols[gene_id], " (", gene_id, ")")
  }
}
rownames(de_expr_matrix) <- row_labels

# Z-score normalization of expression data (rows)
de_expr_z <- t(scale(t(de_expr_matrix)))

# Create sample annotation
sample_anno <- data.frame(CellLine = v$targets$CellLine)
rownames(sample_anno) <- colnames(de_expr_matrix)

# Define annotation colors
ann_colors <- list(
  CellLine = c(Control = "#1B9E77", Y2 = "#D95F02", YT = "#7570B3")
)

# Create heatmap
png("Top_DE_genes_heatmap.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(
  de_expr_z,
  annotation_col = sample_anno,
  annotation_colors = ann_colors,
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  show_rownames = TRUE,
  show_colnames = TRUE,
  main = "Top Differentially Expressed Genes",
  fontsize_row = 8,
  fontsize_col = 10
)
dev.off()

# 2. Heatmap of fold changes across comparisons
# Create a fold change matrix
fc_matrix <- data.frame(
  Y2vsControl = log2FC_Y2vsControl,
  YTvsControl = log2FC_YTvsControl,
  Y2vsYT = log2FC_Y2vsYT
)
fc_matrix <- fc_matrix[all_top_genes, ]
rownames(fc_matrix) <- row_labels

# Create fold change heatmap
png("Fold_change_heatmap.png", width = 8, height = 12, units = "in", res = 300)
pheatmap(
  fc_matrix,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-3, 3, length.out = 101),  # Limit color scale to +/- 3 log2FC
  show_rownames = TRUE,
  main = "Log2 Fold Changes",
  fontsize_row = 8
)
dev.off()

# 3. Sample correlation heatmap
cor_matrix <- cor(norm_expr)
png("Sample_correlation_heatmap.png", width = 10, height = 12, units = "in", res = 300)
pheatmap(
  cor_matrix,
  annotation_col = sample_anno,
  annotation_row = sample_anno,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  display_numbers = TRUE,
  number_format = "%.2f",
  main = "Sample Correlation"
)
dev.off()

# 4. Expression pattern heatmap for specific gene sets (example: cell cycle genes)
# Pathway analysis without prior knowledge of genes in pathways
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# Convert all gene IDs to Entrez IDs for pathway analysis
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
gene_map <- getBM(
  attributes = c("ensembl_gene_id", "entrezgene_id", "external_gene_name"),
  filters = "ensembl_gene_id",
  values = rownames(norm_expr),
  mart = ensembl
)

# Remove entries without Entrez IDs
gene_map <- gene_map[!is.na(gene_map$entrezgene_id),]

# Function to create ranked gene lists for pathway analysis
prepare_gene_list <- function(fc_values, id_map, lfc_threshold = 1) {
  # Sort genes by fold change
  gene_ids <- names(fc_values)
  fc_ordered <- fc_values[order(fc_values, decreasing = TRUE)]
  
  # Map to Entrez IDs, keeping logFC values
  entrez_map <- id_map[match(names(fc_ordered), id_map$ensembl_gene_id),]
  ranked_list <- fc_ordered[!is.na(entrez_map$entrezgene_id)]
  names(ranked_list) <- entrez_map$entrezgene_id[!is.na(entrez_map$entrezgene_id)]
  
  # Extract significant genes based on fold change threshold
  sig_up <- names(ranked_list)[ranked_list >= lfc_threshold]
  sig_down <- names(ranked_list)[ranked_list <= -lfc_threshold]
  
  return(list(
    ranked_list = ranked_list,
    sig_up = sig_up,
    sig_down = sig_down
  ))
}

# Prepare gene lists for each comparison
Y2vsControl_genes <- prepare_gene_list(log2FC_Y2vsControl, gene_map)
YTvsControl_genes <- prepare_gene_list(log2FC_YTvsControl, gene_map)
Y2vsYT_genes <- prepare_gene_list(log2FC_Y2vsYT, gene_map)

# 1. GO Enrichment Analysis
run_go_analysis <- function(gene_list, ont = "BP", title = "") {
  # Run for up-regulated genes
  ego_up <- enrichGO(
    gene = gene_list$sig_up,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Run for down-regulated genes
  ego_down <- enrichGO(
    gene = gene_list$sig_down,
    OrgDb = org.Hs.eg.db,
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  
  # Save results to CSV
  if (nrow(as.data.frame(ego_up)) > 0) {
    write.csv(as.data.frame(ego_up), paste0(title, "_", ont, "_UP.csv"))
    
    # Create and save dotplot
    png(paste0(title, "_", ont, "_UP_dotplot.png"), width = 10, height = 8, units = "in", res = 300)
    print(dotplot(ego_up, showCategory = 15, title = paste0(title, " ", ont, " Up-regulated")))
    dev.off()
  }
  
  if (nrow(as.data.frame(ego_down)) > 0) {
    write.csv(as.data.frame(ego_down), paste0(title, "_", ont, "_DOWN.csv"))
    
    # Create and save dotplot
    png(paste0(title, "_", ont, "_DOWN_dotplot.png"), width = 10, height = 8, units = "in", res = 300)
    print(dotplot(ego_down, showCategory = 15, title = paste0(title, " ", ont, " Down-regulated")))
    dev.off()
  }
  
  return(list(up = ego_up, down = ego_down))
}

# Run GO BP analysis for each comparison
Y2vsControl_GO <- run_go_analysis(Y2vsControl_genes, ont = "BP", title = "Y2vsControl")
YTvsControl_GO <- run_go_analysis(YTvsControl_genes, ont = "BP", title = "YTvsControl")
Y2vsYT_GO <- run_go_analysis(Y2vsYT_genes, ont = "BP", title = "Y2vsYT")

# 2. KEGG Pathway Analysis
run_kegg_analysis <- function(gene_list, title = "") {
  # Run for up-regulated genes
  kk_up <- enrichKEGG(
    gene = gene_list$sig_up,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  # Run for down-regulated genes
  kk_down <- enrichKEGG(
    gene = gene_list$sig_down,
    organism = "hsa",
    pvalueCutoff = 0.05
  )
  
  # Save results to CSV
  if (nrow(as.data.frame(kk_up)) > 0) {
    write.csv(as.data.frame(kk_up), paste0(title, "_KEGG_UP.csv"))
    
    # Create and save dotplot
    png(paste0(title, "_KEGG_UP_dotplot.png"), width = 10, height = 8, units = "in", res = 300)
    print(dotplot(kk_up, showCategory = 15, title = paste0(title, " KEGG Up-regulated")))
    dev.off()
  }
  
  if (nrow(as.data.frame(kk_down)) > 0) {
    write.csv(as.data.frame(kk_down), paste0(title, "_KEGG_DOWN.csv"))
    
    # Create and save dotplot
    png(paste0(title, "_KEGG_DOWN_dotplot.png"), width = 10, height = 8, units = "in", res = 300)
    print(dotplot(kk_down, showCategory = 15, title = paste0(title, " KEGG Down-regulated")))
    dev.off()
  }
  
  return(list(up = kk_up, down = kk_down))
}

# Run KEGG analysis for each comparison
Y2vsControl_KEGG <- run_kegg_analysis(Y2vsControl_genes, title = "Y2vsControl")
YTvsControl_KEGG <- run_kegg_analysis(YTvsControl_genes, title = "YTvsControl")
Y2vsYT_KEGG <- run_kegg_analysis(Y2vsYT_genes, title = "Y2vsYT")

# 3. GSEA Analysis (Gene Set Enrichment Analysis)
run_gsea <- function(gene_list, title = "") {
  gsea_result <- gseGO(
    geneList = gene_list$ranked_list,
    ont = "BP",
    keyType = "ENTREZID",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.05,
    verbose = FALSE,
    OrgDb = org.Hs.eg.db,
    pAdjustMethod = "BH"
  )
  
  # Save results to CSV
  if (nrow(as.data.frame(gsea_result)) > 0) {
    write.csv(as.data.frame(gsea_result), paste0(title, "_GSEA.csv"))
    
    # Create and save GSEA plots for top 5 pathways
    top_pathways <- head(gsea_result@result$ID, 5)
    for (i in seq_along(top_pathways)) {
      png(paste0(title, "_GSEA_", i, ".png"), width = 10, height = 6, units = "in", res = 300)
      print(gseaplot2(gsea_result, geneSetID = i, title = paste0(title, " GSEA")))
      dev.off()
    }
  }
  
  return(gsea_result)
}

# Run GSEA for each comparison
Y2vsControl_GSEA <- run_gsea(Y2vsControl_genes, title = "Y2vsControl")
YTvsControl_GSEA <- run_gsea(YTvsControl_genes, title = "YTvsControl")
Y2vsYT_GSEA <- run_gsea(Y2vsYT_genes, title = "Y2vsYT")

# 4. Create pathway-specific heatmaps for the top pathways from each analysis
create_pathway_heatmap <- function(enrichment_result, gene_map, expr_data, sample_anno, title) {
  if (nrow(as.data.frame(enrichment_result)) > 0) {
    # Get the top pathway
    top_pathway <- enrichment_result@result$ID[1]
    pathway_name <- enrichment_result@result$Description[1]
    
    # Get genes in the pathway
    pathway_genes <- enrichment_result@result$geneID[1]
    pathway_genes <- unlist(strsplit(pathway_genes, "/"))
    
    # Map back to Ensembl IDs
    ensembl_ids <- gene_map$ensembl_gene_id[gene_map$entrezgene_id %in% pathway_genes]
    
    # Get gene symbols for labeling
    gene_names <- gene_map$external_gene_name[gene_map$ensembl_gene_id %in% ensembl_ids]
    gene_labels <- paste0(gene_names, " (", ensembl_ids, ")")
    names(gene_labels) <- ensembl_ids
    
    # Extract expression data for these genes
    pathway_expr <- expr_data[ensembl_ids, ]
    rownames(pathway_expr) <- gene_labels
    
    # Z-score normalize
    pathway_expr_z <- t(scale(t(pathway_expr)))
    
    # Create heatmap
    png(paste0(title, "_top_pathway_heatmap.png"), width = 10, height = 10, units = "in", res = 300)
    pheatmap(
      pathway_expr_z,
      annotation_col = sample_anno,
      annotation_colors = ann_colors,
      main = paste0(title, " - Top Pathway: ", pathway_name),
      fontsize_row = 8
    )
    dev.off()
  }
}

# Create pathway heatmaps for top GO pathways
create_pathway_heatmap(Y2vsControl_GO$up, gene_map, norm_expr, sample_anno, "Y2vsControl_GO_UP")
create_pathway_heatmap(YTvsControl_GO$up, gene_map, norm_expr, sample_anno, "YTvsControl_GO_UP")
create_pathway_heatmap(Y2vsYT_GO$up, gene_map, norm_expr, sample_anno, "Y2vsYT_GO_UP")