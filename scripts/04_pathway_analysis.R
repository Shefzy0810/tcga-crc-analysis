library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)

dds  <- readRDS("data/dds_deseq.rds")
res1 <- results(dds, contrast = c("condition", "tumour", "normal"),
                alpha = 0.05, lfcThreshold = 1)

res1$ensembl <- sub("\\..*", "", rownames(res1))
universe   <- res1$ensembl
up_genes   <- res1$ensembl[which(res1$padj < 0.05 & res1$log2FoldChange > 0)]
down_genes <- res1$ensembl[which(res1$padj < 0.05 & res1$log2FoldChange < 0)]
length(up_genes)
length(down_genes)
ego_up <- enrichGO(gene = up_genes, universe = universe,
                   OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                   ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                   readable = TRUE)

ego_down <- enrichGO(gene = down_genes, universe = universe,
                     OrgDb = org.Hs.eg.db, keyType = "ENSEMBL",
                     ont = "BP", pAdjustMethod = "BH",
                     pvalueCutoff = 0.05, qvalueCutoff = 0.05,
                     readable = TRUE)

head(as.data.frame(ego_up)[, c("Description", "GeneRatio", "p.adjust")], 10)
head(as.data.frame(ego_down)[, c("Description", "GeneRatio", "p.adjust")], 10)
p_up   <- dotplot(ego_up,   showCategory = 15, title = "GO Biological Process: up in tumour")
p_down <- dotplot(ego_down, showCategory = 15, title = "GO Biological Process: down in tumour")
p_up
ggplot2::ggsave("results/go_up_dotplot.png",   p_up,   width = 8, height = 7)
p_down
ggplot2::ggsave("results/go_down_dotplot.png", p_down, width = 8, height = 7)
write.csv(as.data.frame(ego_up),   "results/go_up.csv")
write.csv(as.data.frame(ego_down), "results/go_down.csv")
up_entrez   <- bitr(up_genes,   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
down_entrez <- bitr(down_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
univ_entrez <- bitr(universe,   fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID

kegg_up   <- enrichKEGG(gene = up_entrez,   universe = univ_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_down <- enrichKEGG(gene = down_entrez, universe = univ_entrez, organism = "hsa", pvalueCutoff = 0.05)
kegg_up   <- setReadable(kegg_up,   OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
kegg_down <- setReadable(kegg_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(as.data.frame(kegg_up)[,   c("Description", "GeneRatio", "p.adjust")], 10)
head(as.data.frame(kegg_down)[, c("Description", "GeneRatio", "p.adjust")], 10)
k_up <- dotplot(kegg_up, showCategory = 15, title = "KEGG pathways: up in tumour")
k_up
ggplot2::ggsave("results/kegg_up_dotplot.png", k_up, width = 8, height = 7)
write.csv(as.data.frame(kegg_up),   "results/kegg_up.csv")
write.csv(as.data.frame(kegg_down), "results/kegg_down.csv")
