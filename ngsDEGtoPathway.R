#Setting up environment
install.packages(c("BiocManager"))
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "pathview", "enrichplot", "DOSE"))

install.packages(c("data.table", "ggplot2"))

library(data.table)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(DOSE)
library(enrichplot)

#1. Input Data Format: Import DEG table
deg <- fread("")

p_threshold <- 0.05
fc_threshold <- 2

#Option 1. GO Enrichment Analysis (GSEA)
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

ggsave("dotplot_enrich_go_gsea.png", dotplot_enrich_go_gsea, device = "png", units = "cm", width = 16, height = 18)

#KEGG enrichment
enrich_kegg_gsea <- gseKEGG(geneList = logfc, organism = "hsu")

#Option 2. Over-representation Analysis (ORA)
#GSEA
# make a gene list
deg_up <- deg[logFC > log2(fc_threshold) & adj.P.Val < p_threshold]$V1
# deg_dn <- deg[logFC < -log2(fc_threshold) & adj.P.Val < p_threshold]$V1

enrich_go_fet_up <- enrichGO(gene = deg_up, OrgDb=org.Hs.eg.db, keyType="SYMBOL", ont="ALL", pvalueCutoff=0.05, pAdjustMethod="BH", qvalueCutoff=0.05)

enrich_go_fet_up_df <- enrich_go_fet_up@result

fwrite(enrich_go_fet_up_df, "enrich_go_fet_up_df.tsv", sep = "\t")

#Visualize Top Enriched GO terms
# show the top 10 terms
dotplot_enrich_go_fet_up <- dotplot(enrich_go_fet_up, showCategory = 10, orderBy="GeneRatio")

#KEGG enrichment
enrich_kegg_fet_up <- enrichKEGG(enrich_go_fet_up, organism = "hsu", keyType = "kegg")

ggsave("dotplot_enrich_go_fet_up.png", dotplot_enrich_go_fet_up, device = "png", units = "cm", width = 16, height = 18)