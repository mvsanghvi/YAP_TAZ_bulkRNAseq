if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))
library(data.table)
library(limma)
library(edgeR)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(pheatmap)
library(stats)
library(sva)
setwd("C:/users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq/Rep analysis/Pairwise")

##DATA CLEANING

#Combine count files into single expression matrix
# Retrieve count file paths
file_names <- list.files(path= ".", pattern = "featureCounts_exon\\.txt$", recursive = F, full.names = T)

# Read all count files into a DGElist
dgls <- readDGE(file_names, columns = c(1, 7), skip = 1)
counts_raw <- as.data.frame(dgls$counts)

# Rename samples
sample_name <- gsub("_featureCounts_exon", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sample_name

# Save the combined count table
fwrite(counts_raw, "counts_raw.tsv", sep = "\t", row.names = T)

#Reading raw count data
count_tbl <- fread("counts_raw.tsv", data.table = F)

#Adding Row Names to the Count Table
rownames(count_tbl) <- count_tbl[[1]]
count_tbl <- count_tbl[, -1]

##FILTERING
#Remove genes with little to no expression
perc_keep <- 0.8
gene_keep <- rowSums(count_tbl > 0) >= ceiling(perc_keep * ncol(count_tbl))
count_tbl_low_rm <- count_tbl[gene_keep, ]

#Creating Meta Data table with info about each sample
meta <- data.frame(SampleID = colnames(count_tbl),
                   CellType = c("WT", "WT", "YAPKO", "YAPKO", "YAP_TAZKO", "YAP_TAZKO"))
rownames(meta) <- meta$SampleID
#Combine count data and sample info into object
dge <- DGEList(counts=count_tbl_low_rm, samples = meta)
#Normalize for difficulties between samples
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot=TRUE)

# Use your existing CellType group variable from meta
group <- factor(meta$CellType)
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# Fit the linear model to normalized counts (voom object)
fit <- lmFit(dge_v, design)

# Set up contrasts for pairwise comparisons
contrast.matrix <- makeContrasts(
  YAPKOvsWT = YAPKO - WT,
  YAP_TAZKOvsWT = YAP_TAZKO - WT,
  YAP_TAZKOvsYAPKO = YAP_TAZKO - YAPKO,
  levels = design
)

# Differential Expression Analysis: Compute DE
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract DE results for each comparison
tT_YAPKOvsWT <- topTable(fit2, coef="YAPKOvsWT", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")
tT_YAP_TAZKOvsWT <- topTable(fit2, coef="YAP_TAZKOvsWT", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")
tT_YAP_TAZKOvsYAPKO <- topTable(fit2, coef="YAP_TAZKOvsYAPKO", adjust.method="BH", number=Inf, p.value=1, lfc=0, sort.by = "p")

# Save results
fwrite(tT_YAPKOvsWT, "DE_YAPKO_vs_WT.tsv", sep="\t", row.names=T)
fwrite(tT_YAP_TAZKOvsWT, "DE_YAP_TAZKO_vs_WT.tsv", sep="\t", row.names=T)
fwrite(tT_YAP_TAZKOvsYAPKO, "DE_YAP_TAZKO_vs_YAPKO.tsv", sep="\t", row.names=T)

#Visualizing and GSEA
BiocManager::install(c("pathview", "enrichplot", "DOSE", "msigdbr"))
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
ggsave("dotplot_enrich_go_gsea_2_WT_YK.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

#### FIX THE PATHWAY LIST AND TRY DIFFERENT DATABASES ####

library(msigdbr)
library(clusterProfiler)
library(dplyr)

# 1. FIRST, FIX THE RANKING METRIC TO HANDLE TIES
deg1_ranked <- deg1[order(-logFC)]
deg1_ranked$rank_metric <- sign(deg1_ranked$logFC) * (-log10(deg1_ranked$P.Value))
deg1_ranked <- deg1_ranked[order(-rank_metric)]

gene_list <- deg1_ranked$rank_metric
names(gene_list) <- deg1_ranked$V1
gene_list <- gene_list[!is.na(gene_list) & !is.infinite(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]

cat("Gene list length:", length(gene_list), "\n")
cat("Ties in ranking:", sum(duplicated(gene_list)), "\n")

# 2. TRY HALLMARK PATHWAYS (50 well-defined biological states/processes)
hallmark_pathways <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_list <- hallmark_pathways %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::distinct()

# Check coverage
hallmark_genes_in_data <- sum(unique(hallmark_list$gene_symbol) %in% names(gene_list))
total_hallmark_genes <- length(unique(hallmark_list$gene_symbol))
cat("Hallmark genes in your data:", hallmark_genes_in_data, "/", total_hallmark_genes, "\n")

# Run GSEA with Hallmark
gsea_hallmark <- GSEA(geneList = gene_list,
                      TERM2GENE = hallmark_list,
                      pvalueCutoff = 0.25,
                      minGSSize = 15,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      verbose = TRUE)

if (!is.null(gsea_hallmark) && nrow(as.data.frame(gsea_hallmark)) > 0) {
  gsea_hallmark_df <- as.data.frame(gsea_hallmark)
  print(gsea_hallmark_df[, c("Description", "NES", "pvalue", "p.adjust")])
  fwrite(gsea_hallmark_df, "GSEA_Hallmark_pathways.tsv", sep = "\t", row.names = FALSE)
  
  # Create heatmaps
  expr_matrix <- dge_v$E
  annotation_col <- data.frame(CellType = meta$CellType, row.names = meta$SampleID)
  
  for (i in 1:min(10, nrow(gsea_hallmark_df))) {
    pathway_name <- gsea_hallmark_df$Description[i]
    pathway_genes <- str_split(gsea_hallmark_df$core_enrichment[i], "/")[[1]]
    
    expr_pathway <- expr_matrix[rownames(expr_matrix) %in% pathway_genes, ]
    
    if(nrow(expr_pathway) > 1) {
      pheatmap(expr_pathway,
               scale = "row",
               annotation_col = annotation_col,
               main = str_wrap(pathway_name, width = 50),
               fontsize_row = 8,
               filename = paste0("heatmap_Hallmark_", i, ".png"),
               width = 10,
               height = max(8, nrow(expr_pathway) * 0.25))
    }
  }
} else {
  cat("No significant Hallmark pathways found.\n")
}

# 3. TRY REACTOME PATHWAYS (biological pathways)
reactome_pathways <- msigdbr(species = "Homo sapiens", 
                             category = "C2", 
                             subcategory = "CP:REACTOME")

reactome_list <- reactome_pathways %>% 
  dplyr::select(gs_name, gene_symbol) %>% 
  dplyr::distinct()

reactome_genes_in_data <- sum(unique(reactome_list$gene_symbol) %in% names(gene_list))
total_reactome_genes <- length(unique(reactome_list$gene_symbol))
cat("Reactome genes in your data:", reactome_genes_in_data, "/", total_reactome_genes, "\n")

gsea_reactome <- GSEA(geneList = gene_list,
                      TERM2GENE = reactome_list,
                      pvalueCutoff = 0.25,
                      minGSSize = 15,
                      maxGSSize = 500,
                      pAdjustMethod = "BH",
                      verbose = TRUE)

if (!is.null(gsea_reactome) && nrow(as.data.frame(gsea_reactome)) > 0) {
  gsea_reactome_df <- as.data.frame(gsea_reactome)
  
  # Filter for signaling pathways
  signaling_reactome <- gsea_reactome_df[grepl("WNT|Hippo|FGF|Signaling|MAPK|TGF|Notch|Hedgehog", 
                                               gsea_reactome_df$Description, 
                                               ignore.case = TRUE), ]
  
  print(head(signaling_reactome[, c("Description", "NES", "pvalue", "p.adjust")], 20))
  fwrite(gsea_reactome_df, "GSEA_Reactome_pathways.tsv", sep = "\t", row.names = FALSE)
  
  # Create heatmaps for signaling pathways
  for (i in 1:min(10, nrow(signaling_reactome))) {
    pathway_name <- signaling_reactome$Description[i]
    pathway_genes <- str_split(signaling_reactome$core_enrichment[i], "/")[[1]]
    
    expr_pathway <- expr_matrix[rownames(expr_matrix) %in% pathway_genes, ]
    
    if(nrow(expr_pathway) > 1) {
      pheatmap(expr_pathway,
               scale = "row",
               annotation_col = annotation_col,
               main = str_wrap(pathway_name, width = 50),
               fontsize_row = 7,
               filename = paste0("heatmap_Reactome_signaling_", i, ".png"),
               width = 10,
               height = max(8, nrow(expr_pathway) * 0.25))
    }
  }
} else {
  cat("No significant Reactome pathways found.\n")
}

#### METHOD 3: CREATE SPECIFIC PATHWAY HEATMAPS WITH CUSTOM SETTINGS ####

# Function to create a single pathway heatmap
create_pathway_heatmap <- function(pathway_data, pathway_name_column = "Description", 
                                   expr_matrix, annotation_col, 
                                   output_prefix = "heatmap_pathway") {
  
  pathway_name <- pathway_data[[pathway_name_column]]
  pathway_nes <- pathway_data$NES
  pathway_padj <- pathway_data$p.adjust
  pathway_genes <- str_split(pathway_data$core_enrichment, "/")[[1]]
  
  expr_pathway <- expr_matrix[rownames(expr_matrix) %in% pathway_genes, ]
  
  if(nrow(expr_pathway) < 2) {
    cat("Warning: Less than 2 genes found for", pathway_name, "\n")
    return(NULL)
  }
  
  # Clean up pathway name for display
  display_name <- gsub("HALLMARK_|REACTOME_", "", pathway_name)
  display_name <- gsub("_", " ", display_name)
  
  # Create title
  title <- paste0(display_name, "\n",
                  "NES = ", round(pathway_nes, 2), 
                  ", p.adj = ", format(pathway_padj, digits = 3),
                  " (n=", nrow(expr_pathway), " genes)")
  
  # Clean filename
  file_name_clean <- gsub("[^A-Za-z0-9_]", "_", pathway_name)
  
  # Create heatmap
  pheatmap(expr_pathway,
           scale = "row",
           annotation_col = annotation_col,
           cluster_cols = FALSE,
           clustering_distance_rows = "euclidean",
           main = str_wrap(title, width = 60),
           fontsize_row = 8,
           fontsize_col = 11,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           border_color = NA,
           filename = paste0(output_prefix, "_", file_name_clean, ".png"),
           width = 10,
           height = max(8, nrow(expr_pathway) * 0.3))
  
  cat("Created:", pathway_name, "(", nrow(expr_pathway), "genes )\n")
  
  return(expr_pathway)
}

# Example: Create specific pathway heatmaps

# WNT Signaling
wnt_data <- gsea_hallmark_df[gsea_hallmark_df$Description == "HALLMARK_WNT_BETA_CATENIN_SIGNALING", ]
create_pathway_heatmap(wnt_data, expr_matrix = expr_matrix, 
                       annotation_col = annotation_col, 
                       output_prefix = "heatmap_WNT_signaling")

# Hippo Signaling
hippo_data <- gsea_reactome_df[gsea_reactome_df$Description == "REACTOME_SIGNALING_BY_HIPPO", ]
create_pathway_heatmap(hippo_data, expr_matrix = expr_matrix, 
                       annotation_col = annotation_col, 
                       output_prefix = "heatmap_Hippo_signaling")

# Hedgehog Signaling
hedgehog_data <- gsea_hallmark_df[gsea_hallmark_df$Description == "HALLMARK_HEDGEHOG_SIGNALING", ]
create_pathway_heatmap(hedgehog_data, expr_matrix = expr_matrix, 
                       annotation_col = annotation_col, 
                       output_prefix = "heatmap_Hedgehog_signaling")

# Notch Signaling
notch_data <- gsea_hallmark_df[gsea_hallmark_df$Description == "HALLMARK_NOTCH_SIGNALING", ]
create_pathway_heatmap(notch_data, expr_matrix = expr_matrix, 
                       annotation_col = annotation_col, 
                       output_prefix = "heatmap_Notch_signaling")
# MYC Targets
myc_data <- gsea_hallmark_df[gsea_hallmark_df$Description == "HALLMARK_MYC_TARGETS_V2", ]
create_pathway_heatmap_ordered(myc_data, 
                               expr_matrix = expr_matrix_ordered, 
                               annotation_col = annotation_col_ordered, 
                               output_prefix = "heatmap_ordered_MYC")

# P53 Pathway
p53_data <- gsea_hallmark_df[gsea_hallmark_df$Description == "HALLMARK_P53_PATHWAY", ]
create_pathway_heatmap_ordered(p53_data, 
                               expr_matrix = expr_matrix_ordered, 
                               annotation_col = annotation_col_ordered, 
                               output_prefix = "heatmap_ordered_P53")

##Over-representation Analysis
# make a gene list
deg_up1 <- deg1[logFC > log2(fc_threshold) & adj.P.Val < p_threshold]$V1
# deg_dn <- deg[logFC < -log2(fc_threshold) & adj.P.Val < p_threshold]$V1
enrich_go_fet_up1 <- enrichGO(gene = deg_up1, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="ALL", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)
enrich_go_fet_up_df1 <- enrich_go_fet_up1@result
fwrite(enrich_go_fet_up_df1, "enrich_go_fet_up_df1.tsv", sep = "\t")
# show the top 10 terms
dotplot_enrich_go_fet_up <- dotplot(enrich_go_fet_up, showCategory = 10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_fet_up.png", dotplot_enrich_go_fet_up, device = "png", units = "cm", width = 16, height = 18)

#KEGG Analysis
#enrich_kegg_gsea <- gseKEGG(geneList = logfc, organism = "hsa")
# Prepare ranked gene list
deg_list <- deg1
gene_list <- deg_list$logFC
names(gene_list) <- deg_list$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA or over-representation analysis for GO Biological Process
gsea_res <- gseGO(geneList = gene_list,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  keyType = "ENTREZID",
                  pvalueCutoff = 0.05,
                  verbose = FALSE)
selected_pathways <- c("GO:0008284", "GO:0045787") # example IDs
pathway_genes <- unique(unlist(geneInCategory(gsea_res, selected_pathways)))
expr_subset <- dge_v$E[rownames(dge_v$E) %in% pathway_genes, ]
pheatmap(expr_subset, scale = "row")

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
ggsave("dotplot_enrich_go_gsea_2_WT_YTK.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

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
ggsave("dotplot_enrich_go_gsea_2_YK_YTK.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)