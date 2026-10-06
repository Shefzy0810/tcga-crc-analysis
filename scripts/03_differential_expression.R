library(DESeq2)
dds <- readRDS("data/dds_filtered.rds")
dds <- DESeq(dds)
saveRDS(dds, "data/dds_deseq.rds")
res <- results(dds, contrast = c("condition", "tumour", "normal"), alpha = 0.05)
summary(res)
res1 <- results(dds, contrast = c("condition", "tumour", "normal"),
                alpha = 0.05, lfcThreshold = 1)
summary(res1)
res1$gene_name <- rowData(dds)$gene_name
res_sig <- subset(res1, padj < 0.05)
res_sig <- res_sig[order(res_sig$padj), ]
head(res_sig[, c("gene_name", "log2FoldChange", "padj")], 10)
write.csv(as.data.frame(res_sig), "results/deg_tumour_vs_normal.csv")
library(EnhancedVolcano)
volcano <- EnhancedVolcano(res1,
                           lab = res1$gene_name,
                           x = "log2FoldChange", y = "padj",
                           pCutoff = 0.05, FCcutoff = 1,
                           title = "Tumour vs normal colon (TCGA-COAD)", subtitle = NULL,
                           ylab = bquote(~-Log[10]~ "adjusted P"))

volcano
ggplot2::ggsave("results/volcano_tumour_vs_normal.png", volcano, width = 9, height = 8)
vsd <- vst(dds, blind = FALSE)
up   <- head(res_sig[res_sig$log2FoldChange > 0, ], 25)
down <- head(res_sig[res_sig$log2FoldChange < 0, ], 25)
top  <- rbind(up, down)
mat  <- assay(vsd)[rownames(top), ]
rownames(mat) <- top$gene_name
ann  <- data.frame(condition = dds$condition, row.names = colnames(dds))
library(pheatmap)
pheatmap(mat, scale = "row", annotation_col = ann,
         show_colnames = FALSE, fontsize_row = 6,
         main = "Top 25 up and 25 down genes, tumour vs normal")
pheatmap(mat, scale = "row", annotation_col = ann,
         show_colnames = FALSE, fontsize_row = 6,
         main = "Top 25 up and 25 down genes, tumour vs normal",
         filename = "results/heatmap_top50_degs.png", width = 8, height = 9)

