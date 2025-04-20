# Load necessary libraries
library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)
library(limma)

# For pathway analysis with no replicates, we can:
# a) Use the variance as a ranking metric
# b) Compare cell lines directly using fold changes between pairs

# Option 1a: Use gene variance for GSEA
# Prepare ranked gene list based on variance
gene_list <- gene_var_df$variance
names(gene_list) <- gene_var_df$gene_id
gene_list <- sort(gene_list, decreasing = TRUE)

# Convert Ensembl IDs to Entrez IDs for pathway analysis
id_mapping <- bitr(names(gene_list), 
                   fromType = "ENSEMBL", 
                   toType = "ENTREZID", 
                   OrgDb = org.Hs.eg.db)
  #Getting Warning of failing to map

# Map the ranked values to Entrez IDs
entrez_list <- gene_list[id_mapping$ENSEMBL]
names(entrez_list) <- id_mapping$ENTREZID

# Run GSEA with GO terms
gsea_result <- gseGO(geneList = entrez_list,
                     ont = "BP",  # Biological Process
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     minGSSize = 10,
                     maxGSSize = 500,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH")

# Visualize results
dotplot(gsea_result, showCategory=15, title="GO Enrichment based on Gene Variance")
ggsave("GSEA_variance_GO.png", width=10, height=8)

# Option 1b: Calculate log fold changes between cell lines
# For example, comparing Y2 vs Control
Y2_col <- which(colnames(v$E) == "Y2")
Control_col <- which(colnames(v$E) == "Control")

if(length(Y2_col) > 0 && length(Control_col) > 0) {
  # Calculate simple log fold change
  logFC_Y2_vs_Control <- v$E[, Y2_col] - v$E[, Control_col]
  
  # Create data frame with results
  Y2_vs_Control_df <- data.frame(
    gene_id = rownames(v$E),
    logFC = logFC_Y2_vs_Control
  )
  
  # Sort by absolute fold change
  Y2_vs_Control_df <- Y2_vs_Control_df[order(-abs(Y2_vs_Control_df$logFC)),]
  fwrite(Y2_vs_Control_df, "Y2_vs_Control_logFC.tsv", sep="\t", row.names=FALSE)
  
  # Create ranked list for GSEA
  fc_list <- Y2_vs_Control_df$logFC
  names(fc_list) <- Y2_vs_Control_df$gene_id
  
  # Convert to Entrez IDs
  fc_id_mapping <- bitr(names(fc_list), 
                        fromType = "ENSEMBL", 
                        toType = "ENTREZID", 
                        OrgDb = org.Hs.eg.db)
  
  # Map the ranked values to Entrez IDs
  fc_entrez_list <- fc_list[fc_id_mapping$ENSEMBL]
  # names(fc_entrez_list) <- fc_id_mapping$ENTREZID
  
  # Run GSEA
  gsea_fc_result <- gseGO(geneList = fc_entrez_list,
                          ont = "BP",
                          OrgDb = org.Hs.eg.db,
                          keyType = "ENTREZID",
                          minGSSize = 10,
                          maxGSSize = 500,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
  
  # Visualize results
  dotplot(gsea_fc_result, showCategory=15, title="GO Enrichment Y2 vs Control")
  ggsave("GSEA_Y2vsControl_GO.png", width=10, height=8)
}

# Do the same for YT vs Control
YT_col <- which(colnames(v$E) == "YT")
if(length(YT_col) > 0 && length(Control_col) > 0) {
  # Calculate log fold change
  logFC_YT_vs_Control <- v$E[, YT_col] - v$E[, Control_col]
  
  # Create data frame
  YT_vs_Control_df <- data.frame(
    gene_id = rownames(v$E),
    logFC = logFC_YT_vs_Control
  )
  
  # Sort and save
  YT_vs_Control_df <- YT_vs_Control_df[order(-abs(YT_vs_Control_df$logFC)),]
  fwrite(YT_vs_Control_df, "YT_vs_Control_logFC.tsv", sep="\t", row.names=FALSE)
  
  # GSEA process as above
  # [similar code as for Y2 vs Control]
}