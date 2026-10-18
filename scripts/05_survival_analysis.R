library(DESeq2)
library(survival)
library(survminer)

dds   <- readRDS("data/dds_deseq.rds")
vsd   <- vst(dds, blind = FALSE)
vsd_t <- vsd[, vsd$condition == "tumour"]
clin  <- as.data.frame(colData(vsd_t))

clin$time   <- ifelse(clin$vital_status == "Dead",
                      clin$days_to_death, clin$days_to_last_follow_up)
clin$status <- ifelse(clin$vital_status == "Dead", 1, 0)

table(clin$vital_status, useNA = "ifany")
summary(clin$time)
sum(duplicated(clin$patient))
clin <- clin[clin$vital_status %in% c("Alive", "Dead") &
               !is.na(clin$time) & clin$time > 0, ]
clin <- clin[!duplicated(clin$patient), ]
nrow(clin)
table(clin$status)
deg   <- read.csv("results/deg_tumour_vs_normal.csv", row.names = 1)
genes <- rownames(deg)[1:100]
expr  <- assay(vsd_t)[genes, rownames(clin)]

cox_res <- t(sapply(genes, function(g) {
  fit <- coxph(Surv(clin$time, clin$status) ~ expr[g, ])
  s   <- summary(fit)$coefficients
  c(HR = s[1, "exp(coef)"], p = s[1, "Pr(>|z|)"])
}))
cox_res <- as.data.frame(cox_res)
cox_res$gene_name <- deg[genes, "gene_name"]
cox_res$padj <- p.adjust(cox_res$p, method = "BH")
cox_res <- cox_res[order(cox_res$p), ]
head(cox_res, 10)
write.csv(cox_res, "results/cox_top100_degs.csv")

id <- rownames(cox_res)[cox_res$gene_name == "TEX11"]
clin$TEX11 <- ifelse(expr[id, ] > median(expr[id, ]), "High", "Low")
fit <- survfit(Surv(time, status) ~ TEX11, data = clin)
km <- ggsurvplot(fit, data = clin, pval = TRUE, risk.table = TRUE,
                 title = "TEX11 expression and overall survival")
km

png("results/km_TEX11.png", width = 1800, height = 1600, res = 200)
print(km)
dev.off()
