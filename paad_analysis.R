# Set the working directory
setwd("/Users/Derleen/Desktop/Hackbio_courses/Projects/PAAD_project")
# Set the library path
.libPaths("/Users/Derleen/Desktop/Hackbio_courses/Hackbio_test")
# Install necessary packages 
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("TCGAbiolinks", "edgeR", "limma", "EDASeq", "SummarizedExperiment", "biomaRt"))
install.packages('gplots')
# Load libraries
library(TCGAbiolinks)
library(edgeR)
library(limma)
library(SummarizedExperiment)
library(gplots)
# Step 1: Data Acquisition
# Query TCGA for PAAD RNA-Seq data
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)
# Download data
GDCdownload(query, method = "api", directory = "data")
# Prepare the data
data <- GDCprepare(query, save = TRUE, save.filename = "paad_data.rda", directory = "data")
# Extract expression data and metadata
expression_data <- assay(data)
sample_metadata <- colData(data)
# Step 2: Data Preprocessing
# Convert to DGEList for edgeR
dge <- DGEList(counts = expression_data)

# Create the 'group' variable based on sample_type
dge$samples$group <- factor(ifelse(sample_metadata$sample_type == "Primary Tumor", "Early_Stage", "Normal"))

# Normalize the data
dge <- calcNormFactors(dge)
filtered_dge <- dge[filterByExpr(dge), ]
# Transform the data using voom
v <- voom(filtered_dge, plot = FALSE)
# Step 3: Differential Expression Analysis
# Define the design matrix
design <- model.matrix(~ group, data = dge$samples)
# Fit the model
fit <- lmFit(v, design)
contrast_matrix <- makeContrasts(
  Early_Stage_vs_Normal = 1, -1,
  levels = design
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
# Get differential expression results
results <- topTable(fit2, adjust = "fdr", number = Inf)

# Sort the results by adjusted p-value
results <- results[order(results$adj.P.Val), ]

# Check the top genes by adjusted p-value
head(results)

# Identify significant genes with adjusted p-value < 0.1 and absolute log fold change > 0.5
significant_genes <- results[results$adj.P.Val < 0.1 & abs(results$Early_Stage_vs_Normal) > 0.5, ]

# View the significant genes
head(significant_genes)

# Step 4: Clinical Relevance
# Visualize the results
# Volcano plot
plot(
  x = results$Early_Stage_vs_Normal,
  y = -log10(results$adj.P.Val),
  main = "Volcano Plot for PAAD",
  xlab = "Log Fold Change",
  ylab = "-log10(Adjusted P-Value)",
  pch = 20
)
abline(h = -log10(0.1), col = "red")
abline(v = c(-0.5, 0.5), col = "blue")
# Heatmap of top significant genes
top_genes <- rownames(significant_genes)[1:50]  # Top 50 genes
heatmap.2(
  as.matrix(expression_data[top_genes, ]),
  scale = "row",
  trace = "none",
  dendrogram = "both",
  col = colorRampPalette(c("blue", "white", "red"))(75),
  main = "Heatmap of Top Significant Genes"
)
# Save results for further use
write.csv(significant_genes, "significant_genes.csv")
