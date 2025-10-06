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
ggsave("dotplot_enrich_go_gsea_2_WT_YK.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)
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