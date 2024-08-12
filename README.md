
# Cancer Biology Genomics Project: Pancreatic Adenocarcinoma Analysis

## Overview

This repository contains R scripts for analyzing pancreatic adenocarcinoma (PAAD) RNA-Seq data from The Cancer Genome Atlas (TCGA). The project aims to identify differentially expressed genes between early-stage pancreatic cancer and normal tissue.

## Key Features

- Data acquisition from TCGA
- RNA-Seq data preprocessing
- Differential expression analysis
- Visualization of results (volcano plot and heatmap)
- Export of significant genes for further analysis

## Prerequisites

Ensure you have R installed along with the following packages:

- BiocManager
- TCGAbiolinks
- edgeR
- limma
- EDASeq
- SummarizedExperiment
- biomaRt
- gplots

Install the required packages using:

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("TCGAbiolinks", "edgeR", "limma", "EDASeq", "SummarizedExperiment", "biomaRt"))
install.packages('gplots')

## Usage
1. Run the main script:
  
   paad_analysis.R: https://github.com/Derleen123/Cancer-Biology-Genomics/blob/main/paad_analysis.R
  

## Workflow

1. **Data Acquisition**: Retrieve PAAD RNA-Seq data from TCGA.
2. **Preprocessing**: Normalize data and filter low-expression genes.
3. **Differential Expression Analysis**: Identify significantly differentially expressed genes.
4. **Visualization**: Generate volcano plot and heatmap.
5. **Results Export**: Save significant genes to CSV.

## Output

- `volcano_plot.pdf`: Visualization of differential expression results.
- `heatmap.pdf`: Expression patterns of top significant genes.
- `significant_genes.csv`: List of significantly differentially expressed genes.


## Acknowledgments

This project is part of the Advanced Genomics Course: Bioinformatics for Cancer Biology from HackBio. For more information about the course, visit: https://course.thehackbio.com/
The Cancer Genome Atlas (TCGA) for providing the data. Learn more about TCGA at: https://www.cancer.gov/ccg/research/genome-sequencing/tcga

