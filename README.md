# Tumour vs Normal Gene Expression in Colon Cancer (TCGA-COAD)

## Aim
Identify genes differentially expressed between colon tumours and normal colon tissue, and test whether they relate to patient survival.

## Data
- Source: TCGA-COAD, via the NCI Genomic Data Commons (TCGAbiolinks)
- RNA-seq raw counts (STAR - Counts, unstranded)
- 481 primary tumours, 41 solid tissue normals
- Excluded: 1 metastatic and 1 recurrent sample (too few to analyse)
- Note: groups are unbalanced (481 vs 41), a known TCGA limitation

## Methods
1. **Download and preparation** (`scripts/01_download_data.R`): queried and downloaded STAR gene counts with TCGAbiolinks and merged them into one SummarizedExperiment (60,660 genes × 522 samples).
2. **Filtering and quality control** (`scripts/02_filter_pca.R`): kept genes with at least 10 reads in at least 41 samples (the size of the smallest group), leaving 24,300 genes. Counts were transformed with DESeq2's variance stabilising transformation (VST, blind to group), and the samples were checked by PCA on the 500 most variable genes.
3. **Differential expression** (`scripts/03_differential_expression.R`): DESeq2 with design `~ condition` (normal as reference). Genes were tested against a log2 fold change threshold of 1 (`lfcThreshold = 1`), with significance at Benjamini–Hochberg adjusted p < 0.05. Results were visualised with a volcano plot (EnhancedVolcano) and a heatmap of the top 25 up- and 25 down-regulated genes (pheatmap, row z-scores of VST counts).
4. **Pathway analysis** (`scripts/04_pathway_analysis.R`): over-representation analysis with clusterProfiler, run separately for up- and down-regulated genes, using all 24,300 tested genes as the background (universe). GO Biological Process terms were tested with `enrichGO` (Ensembl IDs) and KEGG pathways with `enrichKEGG` (Ensembl IDs converted to Entrez with `bitr`). Significance was set at BH-adjusted p < 0.05.

## Results
### Quality control
Tumour and normal samples separate clearly along PC1 (23% of variance). Normals cluster tightly, while tumours are more heterogeneous, consistent with known colorectal cancer subtypes. No outliers were removed.

![PCA plot](results/pca_tumour_vs_normal.png)

### Differential expression
**4,821 genes** were differentially expressed (|log2FC| > 1, padj < 0.05): **2,868 up-regulated** and **1,953 down-regulated** in tumours. The full list is in `results/deg_tumour_vs_normal.csv`.

- Top up-regulated genes include CDH3, KRT80, ETV4, SIM2, ESM1, FOXQ1, CEMIP, CLDN1 and WNT2, which point to Wnt signalling, invasion and angiogenesis.
- Down-regulated genes include markers of mature colon epithelium (GUCA2B, AQP8, CLCA4, BEST4, OTOP2, CA7), consistent with loss of normal differentiation in tumours.

![Volcano plot](results/volcano_tumour_vs_normal.png)

![Heatmap](results/heatmap_top50_degs.png)

### Pathway analysis
- **Up in tumour:** KEGG pathways included **Wnt signalling**, the central pathway in colorectal cancer, along with **Cell cycle**, **Homologous recombination**, **Hippo signalling**, and **IL-17 and cytokine signalling**, which point to inflammation in the tumour environment. The top GO terms were developmental and morphogenesis programmes, epithelial differentiation and keratinization, and antimicrobial humoral response.
- **Down in tumour:** normal colon functions (**Mineral absorption**, **Bile secretion**), plus muscle-contraction, neuroactive, calcium and cAMP signalling terms. The latter most likely reflect the smooth muscle and enteric nerve tissue present in normal samples rather than changes inside tumour cells.
- Full tables are in `results/go_up.csv`, `results/go_down.csv`, `results/kegg_up.csv` and `results/kegg_down.csv`.

![GO up](results/go_up_dotplot.png)

![GO down](results/go_down_dotplot.png)

![KEGG up](results/kegg_up_dotplot.png)

## Tools
R 4.6.1, TCGAbiolinks, SummarizedExperiment, DESeq2 1.52.0, EnhancedVolcano 1.30.0, pheatmap 1.0.13, clusterProfiler 4.20.2, org.Hs.eg.db

## Status
In progress: Phases 1–5 complete (setup, data, QC, differential expression, pathway analysis). Next: survival analysis.
