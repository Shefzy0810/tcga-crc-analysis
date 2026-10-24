library(GEOquery); library(limma); library(survival); library(survminer)
options(timeout = 1200)
dir.create("data/GEO", showWarnings = FALSE)

# Download the processed data (about 150-250 MB, run once)
gse  <- getGEO("GSE39582", destdir = "data/GEO", GSEMatrix = TRUE)
eset <- gse[[1]]
dim(eset)
annotation(eset)

# Is the expression already log2-normalised?
summary(exprs(eset)[, 1])

# Which clinical columns do we have?
pd <- pData(eset)
grep("os|rfs|stage|tnm|age|sex", colnames(pd), value = TRUE, ignore.case = TRUE)
# Which samples are tumour vs normal?
table(pd$source_name_ch1)

# Build a clean clinical table
geo <- data.frame(
  time   = as.numeric(pd$`os.delay (months):ch1`),
  status = as.numeric(pd$`os.event:ch1`),
  age    = as.numeric(pd$`age.at.diagnosis (year):ch1`),
  sex    = factor(tolower(pd$`Sex:ch1`), levels = c("female", "male")),
  stage  = pd$`tnm.stage:ch1`,
  row.names = rownames(pd))
table(geo$stage, useNA = "ifany")
table(geo$status, useNA = "ifany")
summary(geo$time)

# Find TEX11 on the microarray
fd <- fData(eset)
tex_probes <- rownames(fd)[fd$`Gene Symbol` == "TEX11"]
tex_probes
round(t(apply(exprs(eset)[tex_probes, ], 1, function(v) c(mean = mean(v), sd = sd(v), max = max(v)))), 2)
geo <- geo[!grepl("non tumoral", pd$source_name_ch1), ]
nrow(geo)

hr_table <- function(fit) {
  s <- summary(fit)
  round(cbind(HR = s$conf.int[, 1], lower95 = s$conf.int[, 3],
              upper95 = s$conf.int[, 4], p = s$coefficients[, 5]), 3)
}

# Prepare the variables
geo$TEX11 <- exprs(eset)["221259_s_at", rownames(geo)]
geo$stage <- factor(c("1" = "I", "2" = "II", "3" = "III", "4" = "IV")[geo$stage],
                    levels = c("I", "II", "III", "IV"))
geo <- geo[!is.na(geo$time) & !is.na(geo$status) & geo$time > 0, ]
nrow(geo); table(geo$status)
saveRDS(geo, "data/geo_clin.rds")

# 1. TEX11 alone
fit_g1 <- coxph(Surv(time, status) ~ TEX11, data = geo)
hr_table(fit_g1)

# 2. TEX11 adjusted for stage, age and sex (same model as TCGA)
fit_g2 <- coxph(Surv(time, status) ~ TEX11 + stage + age + sex, data = geo)
hr_table(fit_g2)
cox.zph(fit_g2)

# TEX11 KM in GEO + save
geo$TEX11_group <- factor(ifelse(geo$TEX11 > median(geo$TEX11), "High", "Low"), levels = c("Low", "High"))
km_g <- ggsurvplot(survfit(Surv(time, status) ~ TEX11_group, data = geo), data = geo,
                   pval = TRUE, risk.table = TRUE, xlab = "Months",
                   title = "TEX11 in GSE39582 (independent validation)")
km_g
png("results/km_TEX11_GSE39582.png", width = 1800, height = 1600, res = 200); print(km_g); dev.off()
write.csv(hr_table(fit_g2), "results/cox_TEX11_GSE39582.csv")

# The 17-gene LASSO signature, applied unchanged
library(SummarizedExperiment)
vsd_t <- readRDS("data/vsd_tumour.rds")
sig   <- readRDS("data/lasso_sig.rds")
w  <- sig$coef * apply(assay(vsd_t)[sig$gene_id, ], 1, sd)
ex <- exprs(eset)[, rownames(geo)]
best_probe <- sapply(sig$gene_name, function(g) {
  p <- rownames(fd)[fd$`Gene Symbol` == g]
  if (length(p) == 0) return(NA)
  p[which.max(rowMeans(ex[p, , drop = FALSE]))]
})
best_probe
found <- !is.na(best_probe)
z <- t(scale(t(ex[best_probe[found], ])))
geo$sig_risk <- as.numeric(t(z) %*% w[found])
hr_table(coxph(Surv(time, status) ~ sig_risk, data = geo))
hr_table(coxph(Surv(time, status) ~ sig_risk + stage + age + sex, data = geo))
write.csv(hr_table(coxph(Surv(time, status) ~ sig_risk + stage + age + sex, data = geo)),
          "results/cox_signature_GSE39582.csv")
write.csv(data.frame(gene = names(best_probe), probe = best_probe),
          "results/signature_probes_GSE39582.csv", row.names = FALSE)
