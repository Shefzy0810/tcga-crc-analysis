library(DESeq2); library(survival); library(survminer); library(glmnet)
vsd_t <- readRDS("data/vsd_tumour.rds")
clin  <- readRDS("data/clin_tumour.rds")
deg   <- read.csv("results/deg_tumour_vs_normal.csv", row.names = 1)

# Candidate genes: protein-coding DEGs (so we can find them on GEO microarrays later)
pc    <- rownames(vsd_t)[rowData(vsd_t)$gene_type == "protein_coding"]
genes <- intersect(rownames(deg), pc)
length(genes)

# x = expression (patients x genes), y = survival
x <- t(assay(vsd_t)[genes, rownames(clin)])
y <- Surv(clin$time, clin$status)

# LASSO-Cox with 10-fold cross-validation
set.seed(2026)
cvfit <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)
plot(cvfit)

# The genes LASSO kept, and their weights
coefs <- coef(cvfit, s = "lambda.min")
keep  <- coefs[, 1] != 0
sig   <- data.frame(gene_id = rownames(coefs)[keep], coef = coefs[keep, 1])
sig$gene_name <- deg[sig$gene_id, "gene_name"]
sig[order(-abs(sig$coef)), ]
nrow(sig)
print(cvfit)

# Risk score = weighted sum of the signature genes
clin$risk <- as.numeric(x[, sig$gene_id] %*% sig$coef)
clin$risk_group <- factor(ifelse(clin$risk > median(clin$risk), "High", "Low"), levels = c("Low", "High"))
fit_km <- survfit(Surv(time, status) ~ risk_group, data = clin)
km <- ggsurvplot(fit_km, data = clin, pval = TRUE, risk.table = TRUE,
                 title = "LASSO signature in TCGA (training set)")
km

write.csv(sig, "results/lasso_signature_genes.csv", row.names = FALSE)
saveRDS(sig, "data/lasso_sig.rds")
# Sensitivity: all expressed protein-coding genes as candidates
x2 <- t(assay(vsd_t)[pc, rownames(clin)])
dim(x2)
set.seed(2026)
cvfit2 <- cv.glmnet(x2, y, family = "cox", alpha = 1, nfolds = 10)
print(cvfit2)
png("results/lasso_cv_degs.png", width = 1600, height = 1200, res = 200); plot(cvfit); dev.off()
png("results/km_lasso_training.png", width = 1800, height = 1600, res = 200); print(km); dev.off()
write.csv(data.frame(candidates = c("DEGs (2,970)", "All protein-coding (16,455)"),
                     cv_deviance_best = c(2.578, 2.601), cv_deviance_null = c(2.598, 2.601),
                     genes_kept = c(17, 10)),
          "results/lasso_cv_summary.csv", row.names = FALSE)
