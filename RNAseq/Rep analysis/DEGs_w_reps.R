if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))

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
                   CellType <- c(rep("WT", 2), rep("YAPKO", 2), rep("YAP_TAZKO", 2)))
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
#plot_save_name <- "PCA_Plot.jpeg"

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
