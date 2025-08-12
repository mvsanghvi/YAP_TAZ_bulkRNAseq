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