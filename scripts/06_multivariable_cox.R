library(DESeq2); library(survival); library(survminer)

# Rebuild the survival table from Chapter 1
dds <- readRDS("data/dds_deseq.rds")
vsd <- vst(dds, blind = FALSE)
vsd_t <- vsd[, vsd$condition == "tumour"]
saveRDS(vsd_t, "data/vsd_tumour.rds")

clin <- as.data.frame(colData(vsd_t))
clin$time <- ifelse(clin$vital_status == "Dead", clin$days_to_death, clin$days_to_last_follow_up)
clin$status <- ifelse(clin$vital_status == "Dead", 1, 0)
clin <- clin[clin$vital_status %in% c("Alive", "Dead") & !is.na(clin$time) & clin$time > 0, ]
clin <- clin[!duplicated(clin$patient), ]
nrow(clin)

# Look at the three new variables
table(clin$ajcc_pathologic_stage, useNA = "ifany")
summary(clin$age_at_index)
table(clin$gender, useNA = "ifany")
# Find the sex column
grep("sex|gender", colnames(clin), value = TRUE, ignore.case = TRUE)
# Simplify stage to I-IV
clin$stage <- gsub("Stage |[ABC]$", "", clin$ajcc_pathologic_stage)
clin$stage <- factor(clin$stage, levels = c("I", "II", "III", "IV"))
table(clin$stage, useNA = "ifany")
# Sex and age as clean variables
table(clin$sex_at_birth, useNA = "ifany")
clin$sex <- factor(tolower(clin$sex_at_birth), levels = c("female", "male"))
clin$age <- clin$age_at_index

# Patients with complete data for the multivariable model
cc <- clin[complete.cases(clin[, c("time", "status", "stage", "age", "sex")]), ]
nrow(cc)
table(cc$status)

saveRDS(clin, "data/clin_tumour.rds")
# 1. Clinical-only model: how much do stage, age and sex matter?
fit_clin <- coxph(Surv(time, status) ~ stage + age + sex, data = cc)
summary(fit_clin)

# 2. Does TEX11 still matter after adjusting for them?
cc$TEX11 <- assay(vsd_t)["ENSG00000120498.14", rownames(cc)]
fit_tex <- coxph(Surv(time, status) ~ TEX11 + stage + age + sex, data = cc)
summary(fit_tex)

# 3. Check the proportional hazards assumption
cox.zph(fit_tex)
anova(fit_clin, fit_tex)
summary(fit_clin)$concordance
hr_table <- function(fit) {
  s <- summary(fit)
  round(cbind(HR = s$conf.int[, 1], lower95 = s$conf.int[, 3],
              upper95 = s$conf.int[, 4], p = s$coefficients[, 5]), 3)
}
write.csv(hr_table(fit_tex), "results/cox_multivariable_TEX11.csv")
fp <- ggforest(fit_tex, data = cc)
ggsave("results/forest_TEX11_multivariable.png", fp, width = 8, height = 5)
fp
