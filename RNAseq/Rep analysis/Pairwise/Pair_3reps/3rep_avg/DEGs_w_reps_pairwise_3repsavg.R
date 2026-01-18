if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "pathview", "enrichplot", "DOSE", "msigdbr"))
install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "tidyverse"))

library(data.table)
library(limma)
library(edgeR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)
library(pheatmap)
library(stringr)
library(dplyr)
library(ggplot2)

setwd("C:/users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq/Rep analysis/Pairwise/Pair_3reps/3rep_avg")

# Data Cleaning and Preprocessing

# 1. Read count files and combine into matrix
file_names <- list.files(path= ".", pattern = "featureCounts_exon\\.txt$", full.names = TRUE)
dgls <- readDGE(file_names, columns = c(1, 7), skip = 1)
counts_raw <- as.data.frame(dgls$counts)

# 2. Rename sample columns
sample_names <- gsub("_featureCounts_exon", "", colnames(counts_raw))
colnames(counts_raw) <- sample_names

# 3. Save combined counts (optional)
fwrite(counts_raw, "counts_raw.tsv", sep = "\t", row.names = TRUE)

# 4. Read raw counts and set rownames
count_tbl <- fread("counts_raw.tsv", data.table = FALSE)
rownames(count_tbl) <- count_tbl[[1]]
count_tbl <- count_tbl[, -1]

# 5. Filter low-expression genes
perc_keep <- 0.8
gene_keep <- rowSums(count_tbl > 0) >= ceiling(perc_keep * ncol(count_tbl))
count_filtered <- count_tbl[gene_keep, ]

# 6. Prepare metadata
meta <- data.frame(SampleID = colnames(count_filtered),
                   CellType = c("WT", "WT", "WT", "YAPKO", "YAPKO", "YAPKO", 
                                "YAP_TAZKO", "YAP_TAZKO", "YAP_TAZKO"))
rownames(meta) <- meta$SampleID

# 7. Create DGEList and normalize
dge <- DGEList(counts = count_filtered, samples = meta)
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot = TRUE)

#Differential Expression Analysis
group <- factor(meta$CellType)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

fit <- lmFit(dge_v, design)

contrast.matrix <- makeContrasts(
  YAPKOvsWT = YAPKO - WT,
  YAP_TAZKOvsWT = YAP_TAZKO - WT,
  YAP_TAZKOvsYAPKO = YAP_TAZKO - YAPKO,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

tT_YAPKOvsWT <- topTable(fit2, coef="YAPKOvsWT", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")
tT_YAP_TAZKOvsWT <- topTable(fit2, coef="YAP_TAZKOvsWT", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")
tT_YAP_TAZKOvsYAPKO <- topTable(fit2, coef="YAP_TAZKOvsYAPKO", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")

fwrite(tT_YAPKOvsWT, "DE_YAPKO_vs_WT.tsv", sep = "\t", row.names = TRUE)
fwrite(tT_YAP_TAZKOvsWT, "DE_YAP_TAZKO_vs_WT.tsv", sep = "\t", row.names = TRUE)
fwrite(tT_YAP_TAZKOvsYAPKO, "DE_YAP_TAZKO_vs_YAPKO.tsv", sep = "\t", row.names = TRUE)
# # Averaging Replicates
# expr_matrix <- dge_v$E   # Normalized expression matrix (log2 CPM)
# meta <- meta[colnames(expr_matrix), ]  # Ensure matching order
# 
# group_factor <- factor(meta$CellType)
# expr_avg <- sapply(levels(group_factor), function(g) {
#   rowMeans(expr_matrix[, group_factor == g, drop = FALSE])
# })
# expr_avg <- as.data.frame(expr_avg)
# rownames(expr_avg) <- rownames(expr_matrix)
# 
# annotation_col_avg <- data.frame(CellType = colnames(expr_avg))
# rownames(annotation_col_avg) <- colnames(expr_avg)

# # Desired order of groups
# desired_order <- c("WT", "YAPKO", "YAP_TAZKO")
# 
# # Reorder columns of averaged expression matrix
# expr_avg <- expr_avg[, desired_order]
# 
# # Reorder annotation to match
# annotation_col_avg <- annotation_col_avg[desired_order,, drop=FALSE]
# 
# # Create ranking metric (e.g., variance across groups) and gene_list
# gene_list <- apply(expr_avg, 1, var)
# names(gene_list) <- rownames(expr_avg)
# gene_list <- gene_list[!is.na(gene_list) & !duplicated(names(gene_list))]
# gene_list <- sort(gene_list, decreasing = TRUE)
# 
# # Get Reactome pathways from msigdbr (subcategory: CP:REACTOME)
# reactome_pathways <- msigdbr(species = "Homo sapiens", 
#                              category = "C2", 
#                              subcategory = "CP:REACTOME") %>%
#                     dplyr::select(gs_name, gene_symbol) %>%
#                     dplyr::distinct()
# 
# # Filter for Notch and Hippo signaling pathways by keywords (case-insensitive)
# notch_pathways <- reactome_pathways %>%
#                   filter(stringr::str_detect(gs_name, regex("Notch", ignore_case = TRUE)))
# 
# hippo_pathways <- reactome_pathways %>%
#                   filter(stringr::str_detect(gs_name, regex("Hippo", ignore_case = TRUE)))
# fgf_pathways <- reactome_pathways %>%
#                   filter(stringr::str_detect(gs_name, regex("FGF", ignore_case = TRUE)))
# 
# gsea_reactome <- GSEA(geneList = gene_list,
#                       TERM2GENE = reactome_pathways,
#                       pvalueCutoff = 0.5,
#                       minGSSize = 15,
#                       maxGSSize = 500,
#                       pAdjustMethod = "BH",
#                       verbose = TRUE)
# 
# gsea_reactome_df <- as.data.frame(gsea_reactome)
# 
# # Filter GSEA results for Notch and Hippo pathways
# notch_gsea <- gsea_reactome_df[stringr::str_detect(gsea_reactome_df$Description, regex("Notch", ignore_case = TRUE)), ]
# hippo_gsea <- gsea_reactome_df[stringr::str_detect(gsea_reactome_df$Description, regex("Hippo", ignore_case = TRUE)), ]
# fgf_gsea <- gsea_reactome_df[stringr::str_detect(gsea_reactome_df$Description, regex("FGF", ignore_case = TRUE)), ]
# # Function create_pathway_heatmap is reused from your pipeline
# 
# # Create heatmaps for Notch pathways
# for(i in seq_len(nrow(notch_gsea))) {
#     create_pathway_heatmap(notch_gsea[i, ], expr_avg, annotation_col_avg, output_prefix = "heatmap_Reactome_Notch")
# }
# 
# # Create heatmaps for Hippo pathways
# for(i in seq_len(nrow(hippo_gsea))) {
#     create_pathway_heatmap(hippo_gsea[i, ], expr_avg, annotation_col_avg, output_prefix = "heatmap_Reactome_Hippo")
# }
# 
# # Create heatmaps for FGF pathways
# for(i in seq_len(nrow(fgf_gsea))) {
#   create_pathway_heatmap(fgf_gsea[i, ], expr_avg, annotation_col_avg, output_prefix = "heatmap_Reactome_FGF")
# }

# # Save averaged matrix for record
# fwrite(data.frame(Gene = rownames(expr_avg), expr_avg), 
#        "expression_matrix_averaged.tsv", sep = "\t", row.names = FALSE)

#Visualizing and GSEA
BiocManager::install(c("pathview", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
deg1 <- fread("DE_YAPKO_vs_WT.tsv")
deg2 <- fread("DE_YAP_TAZKO_vs_WT.tsv")
deg3 <- fread("DE_YAP_TAZKO_vs_YAPKO.tsv")
p_threshold <- 0.05
fc_threshold <- 2

#### Gene Set Enrichment Analysis (GSEA)
### rank the DEGs by the fold change
##WT vs YAPKO
deg1_order_fc <- deg1[order(-logFC)] # rank the genes by logFC in descending order
logfc1 <- deg1_order_fc$logFC # get logFC
names(logfc1) <- deg1_order_fc$V1 # make a named vector of logFC
# GO Analysis
enrich_go_gsea_WT_YK <- gseGO(geneList = logfc1, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, keyType ="SYMBOL", verbose= FALSE)
# save the enrichment table
enrich_go_gsea_df_WT_YK <- enrich_go_gsea_WT_YK@result
fwrite(enrich_go_gsea_df_WT_YK, "enrich_go_gsea_df_WT_YK.tsv", sep = "\t")
#Visualize Top Enriched GO terms
dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea_WT_YK, showCategory = 10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_gsea_2_WT_YKv2.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)


##WT vs YAP/TAZKO
deg2_order_fc <- deg2[order(-logFC)] # rank the genes by logFC in descending order
logfc2 <- deg2_order_fc$logFC # get logFC
names(logfc2) <- deg2_order_fc$V1 # make a named vector of logFC
# GO Analysis
enrich_go_gsea_WT_YTK <- gseGO(geneList = logfc2, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, keyType ="SYMBOL", verbose= FALSE)
# save the enrichment table
enrich_go_gsea_df_WT_YTK <- enrich_go_gsea_WT_YTK@result
fwrite(enrich_go_gsea_df_WT_YTK, "enrich_go_gsea_df_WT_YTK.tsv", sep = "\t")
#Visualize Top Enriched GO terms
dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea_WT_YTK, showCategory = 10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_gsea_2_WT_YTKv2.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

##YAPKO vs YAP/TAZKO
deg3_order_fc <- deg3[order(-logFC)] # rank the genes by logFC in descending order
logfc3 <- deg3_order_fc$logFC # get logFC
names(logfc3) <- deg3_order_fc$V1 # make a named vector of logFC
# GO Analysis
enrich_go_gsea_YK_YTK <- gseGO(geneList = logfc3, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, keyType ="SYMBOL", verbose= FALSE)
# save the enrichment table
enrich_go_gsea_df_YK_YTK <- enrich_go_gsea_YK_YTK@result
fwrite(enrich_go_gsea_df_YK_YTK, "enrich_go_gsea_df_YK_YTK.tsv", sep = "\t")
#Visualize Top Enriched GO terms
dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea_YK_YTK, showCategory = 10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_gsea_2_YK_YTv2K.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

# Prepare Differential Expression for GSEA
# Read differential expression results (example for YAPKO vs WT)
deg <- fread("DE_YAPKO_vs_WT.tsv")
deg$rank_metric <- sign(deg$logFC) * (-log10(deg$P.Value))
deg_ranked <- deg[order(-deg$rank_metric), ]

gene_list <- deg_ranked$rank_metric
names(gene_list) <- deg_ranked$V1
# Remove NAs and duplicates in gene names
valid_idx <- !is.na(gene_list) & !is.infinite(gene_list)
gene_list <- gene_list[valid_idx]
gene_list <- gene_list[!duplicated(names(gene_list))]

#Get Hallmark Pathways from msigdbr 
hallmark_pathways <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::distinct()

# Run GSEA (example)
gsea_hallmark <- GSEA(geneList = gene_list,
                      TERM2GENE = hallmark_pathways,
                      pvalueCutoff = 0.25,
                      minGSSize = 15,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      verbose = TRUE)

gsea_hallmark_df <- as.data.frame(gsea_hallmark)

#Function to create heatmaps with averaged replicates
create_pathway_heatmap <- function(pathway_data, expr_matrix, annotation_col, output_prefix = "heatmap") {
  
  pathway_name <- pathway_data$Description
  pathway_genes <- str_split(pathway_data$core_enrichment, "/")[[1]]
  
  expr_pathway <- expr_matrix[rownames(expr_matrix) %in% pathway_genes, ]
  
  if (nrow(expr_pathway) < 2) {
    message("Warning: Less than 2 genes found for ", pathway_name)
    return(NULL)
  }
  
  file_safe_name <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
  title <- paste0(gsub("HALLMARK_|REACTOME_", "", pathway_name), "\n",
                  "NES = ", round(pathway_data$NES, 2), 
                  ", p.adj = ", format(pathway_data$p.adjust, digits = 3),
                  " (n=", nrow(expr_pathway), " genes)")
  
  pheatmap(expr_pathway,
           scale = "row",
           annotation_col = annotation_col,
           cluster_cols = FALSE,
           fontsize_row = 8,
           main = strwrap(title, width = 60) %>% paste(collapse = "\n"),
           color = colorRampPalette(c("blue", "white", "red"))(100),
           border_color = NA,
           filename = paste0(output_prefix, "_", file_safe_name, ".png"),
           width = 10,
           height = max(8, nrow(expr_pathway) * 0.3))
  
  message("Created heatmap for pathway: ", pathway_name)
}

# --- Example: Generate heatmaps for top 10 Hallmark pathways ---

top_n <- min(10, nrow(gsea_hallmark_df))

for (i in seq_len(top_n)) {
  create_pathway_heatmap(gsea_hallmark_df[i, ], expr_avg, annotation_col_avg, output_prefix = "hallmark_pathway")
}
#### CREATE HEATMAP FOR PLURIPOTENCY GENES ####

library(pheatmap)
library(msigdbr)
library(dplyr)

# Setup
expr_matrix <- dge_v$E
annotation_col <- data.frame(CellType = meta$CellType, row.names = meta$SampleID)
core_pluripotency_genes <- c(
  "POU5F1",    # OCT4
  "SOX2",
  "NANOG",
  "KLF4",
  "MYC",
  "LIN28A",
  "ESRRB",
  "DPPA4",
  "DNMT3B",
  "GDF3",
  "TDGF1",
  "ZFP42",     # REX1
  "UTF1",
  "SALL4",
  "PRDM14",
  "DPPA2",
  "DPPA3",
  "TRIM71",
  "FGF4",
  "TERT",
  "FOXD3",
  "NR0B1"      # DAX1
)

# Check which are in your data
core_genes_in_data <- core_pluripotency_genes[core_pluripotency_genes %in% rownames(expr_matrix)]

cat("Core pluripotency genes found in data:\n")
print(core_genes_in_data)

# Extract expression for these genes
expr_core_pluripotency <- expr_matrix[rownames(expr_matrix) %in% core_genes_in_data, ]
# Check your metadata structure
print(meta)

# Average based on CellType from metadata
expr_averaged <- data.frame(
  WT = rowMeans(expr_core_pluripotency[, meta$CellType == "WT"]),
  YAPKO = rowMeans(expr_core_pluripotency[, meta$CellType == "YAPKO"]),
  YAP_TAZKO = rowMeans(expr_core_pluripotency[, meta$CellType == "YAP_TAZKO"])
)

# Verify dimensions
cat("\nAveraged expression matrix dimensions:\n")
print(dim(expr_averaged))
cat("Genes:", nrow(expr_averaged), "\n")
cat("Conditions:", ncol(expr_averaged), "\n")

# Create annotation for the averaged data
annotation_col_averaged <- data.frame(
  CellType = c("WT", "YAPKO", "YAP_TAZKO"),
  row.names = c("WT", "YAPKO", "YAP_TAZKO")
)

# Create heatmap with averaged data
pheatmap(expr_averaged,
         scale = "row",
         annotation_col = annotation_col_averaged,
         cluster_cols = FALSE,  # Keep the order: WT, YAPKO, YAP_TAZKO
         cluster_rows = TRUE,   # Cluster genes to show similar patterns
         main = "Core Pluripotency Markers\n(Averaged Replicates)",
         fontsize_row = 11,
         fontsize_col = 12,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         border_color = "grey60",
         show_rownames = TRUE,
         show_colnames = TRUE,
         cellwidth = 40,
         cellheight = 20,
         filename = "heatmap_pluripotency_averaged.png",
         width = 6,
         height = 10)
