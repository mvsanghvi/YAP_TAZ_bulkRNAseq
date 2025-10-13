packages_bioc <- c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                   "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                   "pathview", "enrichplot", "DOSE", "msigdbr")
packages_cran <- c("data.table", "readxl", "stringr", "ggplot2", 
                   "ggrepel", "ggfortify", "ggprism", "pheatmap", 
                   "VennDiagram", "corrplot", "Hmisc", "tidyverse")

# Install missing packages
installed <- rownames(installed.packages())
BiocManager::install(setdiff(packages_bioc, installed))
install.packages(setdiff(packages_cran, installed))

# Load all libraries
lapply(c(packages_bioc, packages_cran), library, character.only = TRUE)

# --- Working directory and preprocessing ---
setwd("C:/users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq/Rep analysis/Pairwise/Pair_3reps/3rep_avg")

# --- Averaging replicates ---
# Expression matrix from voom
expr_matrix <- dge_v$E

# Create averaged expression matrix by CellType
celltype_groups <- split(colnames(expr_matrix), meta$CellType)
expr_avg <- sapply(celltype_groups, function(cols) rowMeans(expr_matrix[, cols, drop=FALSE]))
expr_avg <- as.data.frame(expr_avg)

# Annotation for averaged groups
annotation_col_avg <- data.frame(CellType = colnames(expr_avg))
rownames(annotation_col_avg) <- colnames(expr_avg)

# Optional: Save averaged matrix for downstream checks
fwrite(expr_avg, "expression_matrix_averaged.tsv", sep = "\t", row.names = TRUE)

# --- Updated heatmap function using averaged data ---
create_pathway_heatmap <- function(pathway_data, expr_matrix, annotation_col, output_prefix) {
  
  pathway_name <- pathway_data$Description
  pathway_genes <- str_split(pathway_data$core_enrichment, "/")[[1]]
  expr_pathway <- expr_matrix[rownames(expr_matrix) %in% pathway_genes, ]
  
  if (nrow(expr_pathway) < 2) return(NULL)
  
  clean_name <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
  title <- paste0(gsub("HALLMARK_|REACTOME_", "", pathway_name), "\n",
                  "NES = ", round(pathway_data$NES, 2),
                  ", p.adj = ", format(pathway_data$p.adjust, digits = 3),
                  " (n=", nrow(expr_pathway), " genes)")
  
  pheatmap(expr_pathway,
           scale = "row",
           annotation_col = annotation_col,
           cluster_cols = FALSE,
           show_colnames = TRUE,
           main = str_wrap(title, width = 60),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           border_color = NA,
           filename = paste0(output_prefix, "_", clean_name, ".png"))
  cat("Saved averaged heatmap for", pathway_name, "\n")
}

# --- Example of averaged heatmap calls ---
wnt <- subset(gsea_hallmark_df, Description == "HALLMARK_WNT_BETA_CATENIN_SIGNALING")
hippo <- subset(gsea_reactome_df, Description == "REACTOME_SIGNALING_BY_HIPPO")

create_pathway_heatmap(wnt, expr_matrix = expr_avg, annotation_col = annotation_col_avg, output_prefix = "avg_WNT")
create_pathway_heatmap(hippo, expr_matrix = expr_avg, annotation_col = annotation_col_avg, output_prefix = "avg_Hippo")
