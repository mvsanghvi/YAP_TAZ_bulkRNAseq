if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomeInfoDbData")

BiocManager::install(c("sva", "edgeR", "limma", "Biobase", "biomaRt", 
                       "clusterProfiler", "EnhancedVolcano", "org.Hs.eg.db", 
                       "org.Mm.eg.db", "org.Rn.eg.db"))

install.packages(c("data.table", "readxl", "stringr", "ggplot2", "ggrepel", 
                   "ggfortify", "ggprism", "pheatmap", "VennDiagram", 
                   "corrplot", "Hmisc", "stats", "tidyverse"))
library(data.table)
library(edgeR)
library(limma)
library(ggplot2)
library(ggrepel)
library(ggfortify)
library(stats)
library(sva)

getwd()

setwd("C:/Users/mvsan/code/YAP_TAZ_bulkRNAseq/RNAseq")

# Retrieve count file paths
file_names <- list.files(path = ".", pattern = "featureCounts_gene\\.txt$", recursive = T, full.names = T)

# Read all count files into a DGElist
dgls <- readDGE(file_names, columns = c(1, 7), skip = 1)
counts_raw <- as.data.frame(dgls$counts)

# Rename samples
sample_name <- gsub("_featureCounts_gene", "", basename(colnames(counts_raw)))
colnames(counts_raw) <- sample_name

# Save the combined count table
fwrite(counts_raw, "counts_raw.tsv", sep = "\t", row.names = TRUE)

#Read raw count data
count_tbl <- fread("counts_raw.tsv", data.table = F)

#Transform the first column into row names and remove it from the table:
rownames(count_tbl) <- count_tbl[[1]]
count_tbl <- count_tbl[, -1]

#Remove genes with low or no expression (keeps genes that are expressed in at least 80% of the samples):
perc_keep <- 0.8
gene_keep <- rowSums(count_tbl > 0) >= ceiling(perc_keep * ncol(count_tbl))
count_tbl_low_rm <- count_tbl[gene_keep, ]

#Create a metadata table with information about each sample:
meta <- data.frame(SampleID = colnames(count_tbl),
                   CellLine = c("Control", "Y2", "YT"))
rownames(meta) <- meta$SampleID

#Combine the count data and sample information into a DGEList object:

dge <- DGEList(counts=count_tbl_low_rm, samples = meta)

#Normalize the data to adjust for technical differences between samples:

dge <- calcNormFactors(dge, method = "TMM")
dge_v <- voom(dge, plot=TRUE)
# save the dge_v object for later use
saveRDS(dge_v, "dge_v.rds")

#Create PCA plot
shape_column <- "CellLine"
color_column <- "CellLine"
label <- TRUE
label_size <- 4
WT_Y2_YT2 <- "PCA_Plot.pdf"

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

ggsave(WT_Y2_YT2, device = "pdf", units = "cm", width = 16, height = 14)
plot(pca_plot)

#Differential Expression Analysis

#Compare YAP KO to WT control
comparison <- "Y2-Control"

design <- model.matrix(~ 0 + dge_v$targets$CellLine)
colnames(design) <- gsub("dge_v\\$targets\\$CellLine", "", colnames(design))

contrast_matrix <- makeContrasts(
  Y2_vs_Control = Y2 - Control,
  YT_vs_Control = YT - Control,
  YT_vs_Y2 = YT - Y2,
  levels = design
)

v <- voom(dge_v, design)
fit <- lmFit(v, design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

deg_tbl <- topTable(fit_2, coef = colnames(contrast_matrix), n = Inf, p.value=1, lfc=0, sort.by = "p")
fwrite(deg_tbl, "DEG_WT_Y2.tsv", sep = "\t", row.names = T)

#Using DESeq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

library("DESeq2")
library("tidyverse")

#Create DESeq2 Object
condition <- c(rep("Control", 1), rep("Y2", 1), rep("YT", 1))
               
my_colData <- as.data.frame(condition)
rownames(my_colData) <- colnames(data)
my_colData

dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = my_colData,
                              design = ~condition)               
