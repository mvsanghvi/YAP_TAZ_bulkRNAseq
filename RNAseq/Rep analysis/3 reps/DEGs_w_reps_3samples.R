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
setwd("C:/users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq/Rep analysis")

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
                   CellType = c("WT", "WT", "WT", "YAPKO", "YAPKO", "YAPKO", "YAP_TAZKO", "YAP_TAZKO", "YAP_TAZKO"))
rownames(meta) <- meta$SampleID
#Combine count data and sample info into object
dge <- DGEList(counts=count_tbl_low_rm, samples = meta)
#Normalize for difficulties between samples
dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot=TRUE)

##VISULATION
# save the dge_v object for later use
saveRDS(dge_v, "dge_v.rds")

#Create PCA Plot
shape_column <- "CellType"
color_column <- "CellType"
label <- TRUE
label_size <- 4
plot_save_name <- "PCA_Plot.pdf"

meta_table <- dge_v$targets
count_table_t <- as.data.frame(t(dge_v$E))
pca_prep <- prcomp(count_table_t, scale. = TRUE)

pca_plot <- autoplot(pca_prep, label, shape = shape_column, data = meta_table, colour = color_column) +
  geom_text_repel(aes(label = rownames(meta_table)), size = label_size) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor.y=element_blank(),
        panel.grid.major.y=element_blank())

ggsave(plot_save_name, device = "pdf", units = "cm", width = 16, height = 14)

# Read differential expression results
deg_tbl <- fread("DEG_P1_lfc0.tsv")
dge_v <- readRDS("dge_v.rds")

# Set significance thresholds
FC_threshold <- 2
P_threshold <- 0.05

# Create color scheme for different gene categories
keyvals.colour_rd <- ifelse(
  (deg_tbl$logFC < -log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold), 'blue',
  ifelse((deg_tbl$logFC > log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold), 'red', 'grey30')
)

# Add informative labels showing number of genes in each category
names(keyvals.colour_rd)[keyvals.colour_rd == 'blue'] <- paste0(
  "Downregulated (n=", 
  sum(deg_tbl$logFC < -log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold), 
  ")"
)

names(keyvals.colour_rd)[keyvals.colour_rd == 'red'] <- paste0(
  "Upregulated (n=", 
  sum(deg_tbl$logFC > log2(FC_threshold) & deg_tbl$adj.P.Val < P_threshold), 
  ")"
)

names(keyvals.colour_rd)[keyvals.colour_rd == 'grey30'] <- 'Non-significant'


# Prepare sample annotations
col_annot_df <- data.frame(
  CellType = dge_v$targets$CellType,
  row.names = rownames(dge_v$targets)
)
col_annot_df <- col_annot_df[order(col_annot_df$CellType), , drop = FALSE]

# Select significant DEGs
deg_sig <- deg_tbl[abs(deg_tbl$logFC) > log2(FC_threshold) & 
                     deg_tbl$adj.P.Val < P_threshold]$V1
expr_deg <- dge_v$E[deg_sig, rownames(col_annot_df), drop = FALSE]

# Create heatmap
heatmap_custom <- pheatmap(
  expr_deg,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_rownames = FALSE,
  show_colnames = TRUE,
  color = colorRampPalette(c("navy", "white", "red"))(100),
  annotation_col = col_annot_df,
  main = "Differential Expression Heatmap"
)

# Save to image
ggsave("WT_Y2_YT2_rnaseq.png", plot=heatmap_custom)

save_pheatmap_pdf <- function(x, filename, width=8, height=8) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_pdf(heatmap_custom, "heatmap_custom.pdf", width=4, height=4)
# 7. Heatmap of Top Variable Genes
# Calculate variance of each gene
gene_vars <- apply(dge_v$E, 1, var)
gene_var_df <- data.frame(gene_id = names(gene_vars), variance = gene_vars)
gene_var_df <- gene_var_df[order(-gene_var_df$variance),]  # Sort by variance
fwrite(gene_var_df, "gene_variance.tsv", sep="\t", row.names=FALSE)

# Select top 100 most variable genes
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:100]
log_counts_top <- dge_v$E[top_genes, ]

# Scale for heatmap (row-wise)
log_counts_top_scaled <- t(scale(t(log_counts_top)))

# Specify the exact sample order you want
sample_order <- c("WT_1", "WT_2", "WT_3", "Y2_1", "Y2_2", "Y2_3", "YT_1", "YT_2", "YT_3")

# Reorder columns and annotation according to your desired sample order
log_counts_top_scaled_ordered <- log_counts_top_scaled[, sample_order]
annotation_col_ordered <- data.frame(CellType = meta$CellType)
rownames(annotation_col_ordered) <- meta$SampleID
annotation_col_ordered <- annotation_col_ordered[sample_order, , drop = FALSE]

# Plot heatmap with columns ordered and clustering disabled for columns
mymap <- pheatmap(log_counts_top_scaled_ordered,
                  cluster_rows = TRUE,
                  cluster_cols = FALSE,
                  treeheight_row = 0,  # Set dendrogram height to 0 (hides it)
                  treeheight_col = 0,  # Set dendrogram height to 0 (hides it)
                  annotation_col = annotation_col_ordered,
                  show_rownames = TRUE,
                  fontsize_row = 5,
                  fontsize_col = 5,
                  main = "Top 100 Most Variable Genes")

# Save to image
ggsave("WT_Y2_YT2_DEG.png", plot=mymap)
##Differential Expression Analysis
comparison <- "WT-YAPKO-YAP_TAZKO"
design <- model.matrix(~ 0 + dge_v$targets[["CellType"]])
colnames(design) <- gsub(".*]]", "", colnames(design))

contrast_matrix <- makeContrasts(contrasts = comparison, levels = design)

fit <- lmFit(dge_v, design)
fit_2 <- contrasts.fit(fit, contrast_matrix)
fit_2 <- eBayes(fit_2)

deg_tbl <- topTable(fit_2, coef = colnames(contrast_matrix), n = Inf, p.value=1, lfc=0, sort.by = "p")
fwrite(deg_tbl, "DEG_P1_lfc0.tsv", sep = "\t", row.names = T)

BiocManager::install(c("pathview", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
deg <- fread("DEG_P1_lfc0.tsv")

p_threshold <- 0.05
fc_threshold <- 2

## Gene Set Enrichment Analysis
# GO enrichment
# rank the DEGs by the fold change
deg_order_fc <- deg[order(-logFC)] # rank the genes by logFC in descending order
logfc <- deg_order_fc$logFC # get logFC
names(logfc) <- deg_order_fc$V1 # make a named vector of logFC

# GSEA
enrich_go_gsea <- gseGO(geneList = logfc, OrgDb = org.Hs.eg.db, ont = "ALL", pvalueCutoff = 0.05, keyType ="SYMBOL", verbose= FALSE)

# save the enrichment table
enrich_go_gsea_df <- enrich_go_gsea@result

fwrite(enrich_go_gsea_df, "enrich_go_gsea_df.tsv", sep = "\t")

#Visualize Top Enriched GO terms
dotplot_enrich_go_gsea <- dotplot(enrich_go_gsea, showCategory = 10, orderBy="GeneRatio")
ggsave("dotplot_enrich_go_gsea_2.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

##Over-representation Analysis
#GO Enrichment
# make a gene list
deg_up <- deg[logFC > log2(fc_threshold) & adj.P.Val < p_threshold]$V1
# deg_dn <- deg[logFC < -log2(fc_threshold) & adj.P.Val < p_threshold]$V1

enrich_go_fet_up <- enrichGO(gene = deg_up, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="ALL", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)

enrich_go_fet_up_df <- enrich_go_fet_up@result

fwrite(enrich_go_fet_up_df, "enrich_go_fet_up_df.tsv", sep = "\t")
# show the top 10 terms
dotplot_enrich_go_fet_up <- dotplot(enrich_go_fet_up, showCategory = 10, orderBy="GeneRatio")

ggsave("dotplot_enrich_go_fet_up.png", dotplot_enrich_go_fet_up, device = "png", units = "cm", width = 16, height = 18)

##MSigDB Enrichment
