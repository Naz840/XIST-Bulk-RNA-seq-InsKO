# XIST-Bulk-RNA-seq-InsKO
Script for generating PCA, DEG, Heatmap

#!/usr/bin/env Rscript
################################################################################
# RNA-seq Analysis - Customized for Your Metadata
# Using Ins_metadata.csv file
################################################################################

# Install and load required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

required_packages <- c("DESeq2", "ggplot2", "pheatmap", "dplyr", "ggrepel", "RColorBrewer")
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    BiocManager::install(pkg)
    library(pkg, character.only = TRUE)
  }
}

# Create output directories
# WD
getwd()
setwd(r"(C:\Users\shaqu\Box\FMJ lab\Nazmul\Projects\XIST project\RNA seq\Analysis\InsXIST)")
WD <- getwd()
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

################################################################################
# 1. LOAD COUNT DATA
################################################################################

cat("Loading count data...\n")

counts_raw <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/XIST project/RNA seq/Data/InsKO/From FileZilla/03.Result_X202SC25094014-Z01-F001_Mus_musculus/Result_X202SC25094014-Z01-F001_Mus_musculus/3.Quant/1.Count/gene_count.csv", 
                       header = TRUE, 
                       row.names = 1)

View(counts_raw)

# Identify sample columns (those ending with _RNA)
sample_cols <- grep("_RNA$", colnames(counts_raw))
count_data <- counts_raw[, sample_cols]
rownames(count_data) <- counts_raw$gene_id

# Extract gene annotations
annotation_cols <- setdiff(1:ncol(counts_raw), c(1, sample_cols))
gene_annotations <- counts_raw[, annotation_cols]
rownames(gene_annotations) <- counts_raw$gene_id

cat("Total genes:", nrow(count_data), "\n")
cat("Total samples in count file:", ncol(count_data), "\n")

################################################################################
# 2. LOAD YOUR METADATA
################################################################################

cat("\nLoading your metadata file...\n")

# Read metadata
metadata <- read.csv("C:/Users/shaqu/Box/FMJ lab/Nazmul/Projects/XIST project/RNA seq/Data/InsKO/From FileZilla/03.Result_X202SC25094014-Z01-F001_Mus_musculus/Result_X202SC25094014-Z01-F001_Mus_musculus/3.Quant/1.Count/Ins_metadata.csv", stringsAsFactors = FALSE)

cat("Metadata loaded with", nrow(metadata), "samples\n")
cat("Columns in metadata:", paste(colnames(metadata), collapse = ", "), "\n")

# Set SampleName as rownames
rownames(metadata) <- metadata$SampleName

# Check if all samples in count data have metadata
if (!all(colnames(count_data) %in% metadata$SampleName)) {
  missing <- setdiff(colnames(count_data), metadata$SampleName)
  cat("\nWARNING: These samples in count data are missing from metadata:\n")
  print(missing)
}

# Check if all metadata samples are in count data
if (!all(metadata$SampleName %in% colnames(count_data))) {
  missing <- setdiff(metadata$SampleName, colnames(count_data))
  cat("\nWARNING: These samples in metadata are missing from count data:\n")
  print(missing)
}

# Keep only samples present in both
common_samples <- intersect(colnames(count_data), metadata$SampleName)
count_data <- count_data[, common_samples]
metadata <- metadata[common_samples, ]

cat("\nFinal number of samples:", length(common_samples), "\n")

# Convert Sex to full names for better plots
metadata$SexFull <- ifelse(metadata$Sex == "F", "Female", "Male")

# Convert to factors with proper ordering
metadata$Genotype <- factor(metadata$Genotype, levels = c("WT", "InsXistKO"))
metadata$SexFull <- factor(metadata$SexFull, levels = c("Female", "Male"))
metadata$Status <- factor(metadata$Status, 
                          levels = c("Female_InsWT", "Female_InsKO", 
                                     "Male_InsWT", "Male_InsKO"))

# Display metadata summary
cat("\n=== Metadata Summary ===\n")
print(head(metadata))
cat("\nSamples per Status:\n")
print(table(metadata$Status))
cat("\nSamples per Genotype:\n")
print(table(metadata$Genotype))
cat("\nSamples per Sex:\n")
print(table(metadata$SexFull))

################################################################################
# 3. PREPARE COUNT MATRIX
################################################################################

cat("\nPreparing count matrix...\n")

# Convert to integer matrix
count_matrix <- as.matrix(count_data)
mode(count_matrix) <- "integer"

# Filter low count genes (at least 10 total reads)
keep <- rowSums(count_matrix) >= 10
count_matrix_filtered <- count_matrix[keep, ]
gene_annotations_filtered <- gene_annotations[keep, ]

cat("Genes before filtering:", nrow(count_matrix), "\n")
cat("Genes after filtering:", nrow(count_matrix_filtered), "\n\n")

################################################################################
# 4. DESeq2 ANALYSIS
################################################################################

cat("Running DESeq2 analysis...\n")

# Create DESeq2 dataset
# Design: test for Genotype effect while accounting for Sex
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix_filtered,
  colData = metadata,
  design = ~ SexFull + Genotype
)

cat("Design formula: ~ SexFull + Genotype\n")
cat("This tests the Genotype effect (InsXistKO vs WT) while controlling for Sex differences\n\n")

# Run DESeq2
dds <- DESeq(dds)

# Get normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "results/normalized_counts.csv")

################################################################################
# 5. DIFFERENTIAL EXPRESSION RESULTS
################################################################################

cat("Extracting differential expression results...\n")

# Main comparison: InsXistKO vs WT (controlling for sex)
res <- results(dds, contrast = c("Genotype", "InsXistKO", "WT"), alpha = 0.05)

# Convert to data frame and add annotations
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)
res_df$gene_name <- gene_annotations_filtered[rownames(res_df), "gene_name"]
res_df$gene_biotype <- gene_annotations_filtered[rownames(res_df), "gene_biotype"]
res_df$gene_description <- gene_annotations_filtered[rownames(res_df), "gene_description"]

# Reorder columns
res_df <- res_df[, c("gene_id", "gene_name", "gene_biotype", "baseMean", 
                     "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", 
                     "gene_description")]

# Sort by adjusted p-value
res_df <- res_df[order(res_df$padj), ]

# Save results
write.csv(res_df, "results/DEG_InsXistKO_vs_WT_all_genes.csv", row.names = FALSE)

# Save only significant genes
res_sig <- res_df[res_df$padj < 0.05 & !is.na(res_df$padj), ]
write.csv(res_sig, "results/DEG_InsXistKO_vs_WT_significant.csv", row.names = FALSE)

# Print summary
cat("\n=== DIFFERENTIAL EXPRESSION SUMMARY ===\n")
cat("Total genes tested:", nrow(res_df), "\n")
cat("Significant genes (padj < 0.05):", sum(res_df$padj < 0.05, na.rm = TRUE), "\n")
cat("  Upregulated in InsXistKO:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE), "\n")
cat("  Downregulated in InsXistKO:", 
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE), "\n")
cat("Strong effect (padj < 0.05 & |log2FC| > 1):", 
    sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE), "\n")
cat("Very strong effect (padj < 0.01 & |log2FC| > 2):", 
    sum(res_df$padj < 0.01 & abs(res_df$log2FoldChange) > 2, na.rm = TRUE), "\n\n")

cat("Top 15 differentially expressed genes:\n")
print(res_df[1:15, c("gene_name", "log2FoldChange", "padj")])
cat("\n")

################################################################################
# 6. SEX-SPECIFIC ANALYSES (OPTIONAL)
################################################################################

cat("Running sex-specific analyses...\n")

# Female-only analysis
dds_female <- dds[, dds$SexFull == "Female"]
dds_female$Genotype <- droplevels(dds_female$Genotype)
res_female <- results(dds_female, contrast = c("Genotype", "InsXistKO", "WT"))
res_female_df <- as.data.frame(res_female)
res_female_df$gene_name <- gene_annotations_filtered[rownames(res_female_df), "gene_name"]
res_female_df <- res_female_df[order(res_female_df$padj), ]
write.csv(res_female_df, "results/DEG_Female_InsXistKO_vs_WT.csv", row.names = TRUE)

cat("Female-only DEGs (padj < 0.05):", sum(res_female_df$padj < 0.05, na.rm = TRUE), "\n")

# Male-only analysis
dds_male <- dds[, dds$SexFull == "Male"]
dds_male$Genotype <- droplevels(dds_male$Genotype)
res_male <- results(dds_male, contrast = c("Genotype", "InsXistKO", "WT"))
res_male_df <- as.data.frame(res_male)
res_male_df$gene_name <- gene_annotations_filtered[rownames(res_male_df), "gene_name"]
res_male_df <- res_male_df[order(res_male_df$padj), ]
write.csv(res_male_df, "results/DEG_Male_InsXistKO_vs_WT.csv", row.names = TRUE)

cat("Male-only DEGs (padj < 0.05):", sum(res_male_df$padj < 0.05, na.rm = TRUE), "\n\n")

################################################################################
# 7. PCA PLOT
################################################################################

cat("Generating PCA plot...\n")

# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# Calculate PCA
pca_data <- plotPCA(vsd, intgroup = c("SexFull", "Genotype"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# PCA plot colored by Genotype, shaped by Sex
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, 
                                 color = Genotype, 
                                 shape = SexFull)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
  scale_color_manual(values = c("WT" = "#3B9AB2", "InsXistKO" = "#E8601C")) +
  scale_shape_manual(values = c("Female" = 16, "Male" = 17)) +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    title = "PCA of RNA-seq Samples",
    color = "Genotype",
    shape = "Sex"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    aspect.ratio = 1
  )

ggsave("figures/PCA_plot.png", pca_plot, width = 8, height = 6, dpi = 300)
ggsave("figures/PCA_plot.pdf", pca_plot, width = 8, height = 6)

# PCA colored by Status (combined groups)
pca_data_status <- plotPCA(vsd, intgroup = "Status", returnData = TRUE)

pca_plot_status <- ggplot(pca_data_status, aes(x = PC1, y = PC2, color = Status)) +
  geom_point(size = 4, alpha = 0.8) +
  geom_text_repel(aes(label = name), size = 3, max.overlaps = 20) +
  scale_color_brewer(palette = "Set2") +
  labs(
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance"),
    title = "PCA by Status Group",
    color = "Status"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    aspect.ratio = 1
  )

ggsave("figures/PCA_plot_by_status.png", pca_plot_status, width = 8, height = 6, dpi = 300)
ggsave("figures/PCA_plot_by_status.pdf", pca_plot_status, width = 8, height = 6)

################################################################################
# 8. VOLCANO PLOT
################################################################################

cat("Generating volcano plot...\n")

res_plot <- res_df[!is.na(res_df$padj), ]

# Add significance categories
res_plot$significance <- "Not Sig"
res_plot$significance[res_plot$padj < 0.05 & res_plot$log2FoldChange > 1] <- "Up"
res_plot$significance[res_plot$padj < 0.05 & res_plot$log2FoldChange < -1] <- "Down"
res_plot$significance[res_plot$padj < 0.05 & 
                        abs(res_plot$log2FoldChange) <= 1] <- "Sig"

# Count genes in each category
cat("\nVolcano plot statistics:\n")
cat("Upregulated (padj<0.05, log2FC>1):", sum(res_plot$significance == "Up"), "\n")
cat("Downregulated (padj<0.05, log2FC<-1):", sum(res_plot$significance == "Down"), "\n")
cat("Significant (padj<0.05, |log2FC|<=1):", sum(res_plot$significance == "Sig"), "\n")

# Label top genes
res_plot$label <- ""
top_up <- head(res_plot[res_plot$significance == "Up", ], 10)
top_down <- head(res_plot[res_plot$significance == "Down", ], 10)

if(nrow(top_up) > 0) {
  res_plot$label[match(top_up$gene_id, res_plot$gene_id)] <- 
    ifelse(is.na(top_up$gene_name) | top_up$gene_name == "", 
           top_up$gene_id, top_up$gene_name)
}
if(nrow(top_down) > 0) {
  res_plot$label[match(top_down$gene_id, res_plot$gene_id)] <- 
    ifelse(is.na(top_down$gene_name) | top_down$gene_name == "", 
           top_down$gene_id, top_down$gene_name)
}

volcano <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = significance), alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c("Up" = "#E8601C", 
               "Down" = "#3B9AB2", 
               "Sig" = "#90C987",
               "Not Sig" = "grey70"),
    labels = c("Up" = "Upregulated", 
               "Down" = "Downregulated", 
               "Sig" = "Significant",
               "Not Sig" = "Not Significant")
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  geom_text_repel(aes(label = label), size = 3, max.overlaps = 20,
                  box.padding = 0.5, point.padding = 0.3,
                  min.segment.length = 0) +
  labs(
    title = "Differential Gene Expression: InsXistKO vs WT",
    x = "log2 Fold Change",
    y = "-log10 (adjusted p-value)",
    color = "Category"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3)
  )

ggsave("figures/volcano_plot.png", volcano, width = 10, height = 8, dpi = 300)
ggsave("figures/volcano_plot.pdf", volcano, width = 10, height = 8)

################################################################################
# 9. HEATMAP OF TOP DEGs
################################################################################

cat("Generating heatmap of top DEGs...\n")

# Get top significant genes
top_n <- min(50, sum(res_df$padj < 0.05, na.rm = TRUE))

if(top_n > 0) {
  top_genes <- head(res_df[!is.na(res_df$padj), ], top_n)
  top_counts <- normalized_counts[rownames(top_genes), ]
  
  # Z-score normalization
  top_counts_scaled <- t(scale(t(top_counts)))
  
  # Use gene names for rownames
  rownames(top_counts_scaled) <- ifelse(
    is.na(top_genes$gene_name) | top_genes$gene_name == "",
    rownames(top_genes),
    top_genes$gene_name
  )
  
  # Create annotation
  annotation_col <- data.frame(
    Genotype = metadata[colnames(top_counts), "Genotype"],
    Sex = metadata[colnames(top_counts), "SexFull"],
    Status = metadata[colnames(top_counts), "Status"],
    row.names = colnames(top_counts)
  )
  
  # Define colors
  ann_colors <- list(
    Genotype = c(WT = "#3B9AB2", InsXistKO = "#E8601C"),
    Sex = c(Female = "#F21A00", Male = "#3B9AB2"),
    Status = c(Female_InsWT = "#78B7C5", 
               Female_InsKO = "#EBCC2A",
               Male_InsWT = "#3B9AB2",
               Male_InsKO = "#E8601C")
  )
  
  # Create heatmap
  png("figures/heatmap_top_DEGs.png", width = 12, height = 10, units = "in", res = 300)
  pheatmap(top_counts_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = paste("Top", top_n, "Differentially Expressed Genes"),
           fontsize_row = 8,
           fontsize_col = 9,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           border_color = NA)
  dev.off()
  
  pdf("figures/heatmap_top_DEGs.pdf", width = 12, height = 10)
  pheatmap(top_counts_scaled,
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           show_rownames = TRUE,
           show_colnames = TRUE,
           main = paste("Top", top_n, "Differentially Expressed Genes"),
           fontsize_row = 8,
           fontsize_col = 9,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           border_color = NA)
  dev.off()
  
  cat("Heatmap created with", top_n, "genes\n")
} else {
  cat("No significant genes found for heatmap\n")
}

################################################################################
# 10. SAMPLE DISTANCE HEATMAP
################################################################################

cat("Generating sample distance heatmap...\n")

sample_dists <- dist(t(assay(vsd)))
sample_dist_matrix <- as.matrix(sample_dists)

# Shorten sample names for display
colnames(sample_dist_matrix) <- gsub("_RNA", "", colnames(sample_dist_matrix))
rownames(sample_dist_matrix) <- gsub("_RNA", "", rownames(sample_dist_matrix))

annotation_col_dist <- data.frame(
  Genotype = metadata$Genotype,
  Sex = metadata$SexFull,
  row.names = gsub("_RNA", "", rownames(metadata))
)

png("figures/sample_distance_heatmap.png", width = 10, height = 9, units = "in", res = 300)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col_dist,
         annotation_row = annotation_col_dist,
         annotation_colors = list(
           Genotype = c(WT = "#3B9AB2", InsXistKO = "#E8601C"),
           Sex = c(Female = "#F21A00", Male = "#3B9AB2")
         ),
         main = "Sample-to-Sample Distance",
         fontsize = 10)
dev.off()

pdf("figures/sample_distance_heatmap.pdf", width = 10, height = 9)
pheatmap(sample_dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col_dist,
         annotation_row = annotation_col_dist,
         annotation_colors = list(
           Genotype = c(WT = "#3B9AB2", InsXistKO = "#E8601C"),
           Sex = c(Female = "#F21A00", Male = "#3B9AB2")
         ),
         main = "Sample-to-Sample Distance",
         fontsize = 10)
dev.off()

################################################################################
# 11. MA PLOT
################################################################################

cat("Generating MA plot...\n")

png("figures/MA_plot.png", width = 8, height = 6, units = "in", res = 300)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: InsXistKO vs WT")
dev.off()

pdf("figures/MA_plot.pdf", width = 8, height = 6)
plotMA(res, ylim = c(-5, 5), main = "MA Plot: InsXistKO vs WT")
dev.off()

################################################################################
# 12. CREATE SUMMARY REPORT
################################################################################

summary_stats <- data.frame(
  Comparison = c("All samples: InsXistKO vs WT",
                 "Female only: InsXistKO vs WT",
                 "Male only: InsXistKO vs WT"),
  Total_Genes_Tested = c(nrow(res_df), nrow(res_female_df), nrow(res_male_df)),
  Significant_padj_0.05 = c(
    sum(res_df$padj < 0.05, na.rm = TRUE),
    sum(res_female_df$padj < 0.05, na.rm = TRUE),
    sum(res_male_df$padj < 0.05, na.rm = TRUE)
  ),
  Upregulated = c(
    sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE),
    sum(res_female_df$padj < 0.05 & res_female_df$log2FoldChange > 0, na.rm = TRUE),
    sum(res_male_df$padj < 0.05 & res_male_df$log2FoldChange > 0, na.rm = TRUE)
  ),
  Downregulated = c(
    sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE),
    sum(res_female_df$padj < 0.05 & res_female_df$log2FoldChange < 0, na.rm = TRUE),
    sum(res_male_df$padj < 0.05 & res_male_df$log2FoldChange < 0, na.rm = TRUE)
  )
)

write.csv(summary_stats, "results/analysis_summary.csv", row.names = FALSE)

cat("\n=== Analysis Summary ===\n")
print(summary_stats)

################################################################################
# COMPLETE
################################################################################

cat("\n================================================\n")
cat("✓ Analysis Complete!\n")
cat("================================================\n")
cat("\nResults folder contains:\n")
cat("  - DEG_InsXistKO_vs_WT_all_genes.csv\n")
cat("  - DEG_InsXistKO_vs_WT_significant.csv\n")
cat("  - DEG_Female_InsXistKO_vs_WT.csv\n")
cat("  - DEG_Male_InsXistKO_vs_WT.csv\n")
cat("  - normalized_counts.csv\n")
cat("  - analysis_summary.csv\n")
cat("\nFigures folder contains:\n")
cat("  - PCA_plot (by genotype and sex)\n")
cat("  - PCA_plot_by_status\n")
cat("  - volcano_plot\n")
cat("  - heatmap_top_DEGs\n")
cat("  - sample_distance_heatmap\n")
cat("  - MA_plot\n")
cat("\nAll figures saved in PNG (300 DPI) and PDF formats\n")
cat("================================================\n")
