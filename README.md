# Tumour vs Normal Gene Expression in Colon Cancer (TCGA-COAD)

## Aim
Identify genes differentially expressed between colon tumours and normal colon tissue, and test whether they relate to patient survival.

## Data
- Source: TCGA-COAD, via the NCI Genomic Data Commons (TCGAbiolinks)
- RNA-seq raw counts (STAR - Counts)
- 481 primary tumours, 41 solid tissue normals
- Excluded: 1 metastatic and 1 recurrent sample (too few to analyse)
- Note: groups are unbalanced (481 vs 41), a known TCGA limitation

## Tools
R 4.6.1, DESeq2 1.52.0, TCGAbiolinks

## Status
In progress
