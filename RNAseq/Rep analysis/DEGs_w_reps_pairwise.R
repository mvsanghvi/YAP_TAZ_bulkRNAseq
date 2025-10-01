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
tt_YAPKOvsWT <- topTable(fit2, coef="YAPKOvsWT", adjust.method="BH", number=Inf)
tt_YAP_TAZKOvsWT <- topTable(fit2, coef="YAP_TAZKOvsWT", adjust.method="BH", number=Inf)
tt_YAP_TAZKOvsYAPKO <- topTable(fit2, coef="YAP_TAZKOvsYAPKO", adjust.method="BH", number=Inf)

# Save results
fwrite(tt_YAPKOvsWT, "DE_YAPKO_vs_WT.tsv", sep="\t", row.names=T)
fwrite(tt_YAP_TAZKOvsWT, "DE_YAP_TAZKO_vs_WT.tsv", sep="\t", row.names=T)
fwrite(tt_YAP_TAZKOvsYAPKO, "DE_YAP_TAZKO_vs_YAPKO.tsv", sep="\t", row.names=T)
