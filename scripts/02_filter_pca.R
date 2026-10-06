library(SummarizedExperiment)
library(DESeq2)
se <- readRDS("data/tcga_coad_se.rds")
table(se$sample_type)
se$condition <- factor(
  ifelse(se$sample_type == "Solid Tissue Normal", "normal", "tumour"),
  levels = c("normal", "tumour")
)
table(se$condition)
assays(se) <- assays(se)["unstranded"]
dds <- DESeqDataSet(se, design = ~ condition)
dds
keep <- rowSums(counts(dds) >= 10) >= 41
dds <- dds[keep, ]
nrow(dds)
vsd <- vst(dds, blind = TRUE)
pca_plot <- plotPCA(vsd, intgroup = "condition")
pca_plot
ggplot2::ggsave("results/pca_tumour_vs_normal.png", pca_plot, width = 7, height = 5)
