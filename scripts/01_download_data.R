# Step 1: Find TCGA colon cancer RNA-seq data
library(TCGAbiolinks)

query <- GDCquery(
  project = "TCGA-COAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

results <- getResults(query)
table(results$sample_type)

# Step 2: Download the files (about 2 GB, run once)
GDCdownload(query, directory = "data/GDCdata", files.per.chunk = 50)
