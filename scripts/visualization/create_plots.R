#!/usr/bin/env Rscript
# =============================================================================
# Visualization Script for NGS Bioinformatics Project
# Course: LIFE622C - Next Generation Sequencing and Bioinformatics
# Institution: POSTECH
# Student: Saggurthi Dinesh (BE21B032)
# Faculty: Dr. Jong Kyoung Kim
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(VennDiagram)
  library(RColorBrewer)
  library(corrplot)
  library(pcaMethods)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
  library(cowplot)
})

# =============================================================================
# Configuration and Setup
# =============================================================================

# Set working directory and create output directories
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(results_dir, "figures")

# Create figures directory if it doesn't exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# Color palette for plots
color_palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# =============================================================================
# Helper Functions
# =============================================================================

# Function to save plots with consistent formatting
save_plot <- function(plot_obj, filename, width = 10, height = 8, dpi = 300) {
  full_path <- file.path(figures_dir, filename)
  
  if (inherits(plot_obj, "ggplot")) {
    ggsave(full_path, plot_obj, width = width, height = height, dpi = dpi)
  } else {
    # For base R plots or other plot types
    png(full_path, width = width * 100, height = height * 100, res = dpi)
    print(plot_obj)
    dev.off()
  }
  
  cat("Saved plot:", full_path, "\n")
}

# Function to create a sample correlation heatmap
create_correlation_heatmap <- function(data_matrix, filename = "sample_correlation_heatmap.png") {
  # Calculate correlation matrix
  cor_matrix <- cor(data_matrix, use = "complete.obs")
  
  # Create heatmap
  pheatmap(
    cor_matrix,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    show_colnames = TRUE,
    fontsize = 8,
    main = "Sample Correlation Matrix",
    filename = file.path(figures_dir, filename),
    width = 10,
    height = 8
  )
  
  cat("Created correlation heatmap:", filename, "\n")
}

# Function to create PCA plot
create_pca_plot <- function(data_matrix, sample_groups = NULL, filename = "pca_plot.png") {
  # Perform PCA
  pca_result <- prcomp(t(data_matrix), scale. = TRUE)
  
  # Prepare data for plotting
  pca_data <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = rownames(pca_result$x)
  )
  
  if (!is.null(sample_groups)) {
    pca_data$Group <- sample_groups
  }
  
  # Create PCA plot
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_text_repel(aes(label = Sample), size = 3, max.overlaps = 10) +
    theme_minimal() +
    labs(
      title = "Principal Component Analysis",
      x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100, 1), "%)"),
      y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100, 1), "%)")
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  if (!is.null(sample_groups)) {
    pca_plot <- pca_plot + aes(color = Group) + scale_color_brewer(palette = "Set1")
  }
  
  save_plot(pca_plot, filename)
  return(pca_plot)
}

# Function to create volcano plot
create_volcano_plot <- function(deg_data, filename = "volcano_plot.png") {
  # Assuming deg_data has columns: gene, log2FoldChange, padj
  # Create significance threshold
  deg_data$significant <- ifelse(deg_data$padj < 0.05 & abs(deg_data$log2FoldChange) > 1, "Significant", "Not Significant")
  
  # Create volcano plot
  volcano_plot <- ggplot(deg_data, aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.6, size = 1) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
    theme_minimal() +
    labs(
      title = "Volcano Plot - Differential Expression",
      x = "log2 Fold Change",
      y = "-log10 Adjusted P-value",
      color = "Significance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom"
    )
  
  save_plot(volcano_plot, filename)
  return(volcano_plot)
}

# Function to create Manhattan plot
create_manhattan_plot <- function(eqtl_data, filename = "manhattan_plot.png") {
  # Assuming eqtl_data has columns: chromosome, position, pvalue
  # Create Manhattan plot
  manhattan_plot <- ggplot(eqtl_data, aes(x = position, y = -log10(pvalue), color = factor(chromosome))) +
    geom_point(alpha = 0.6, size = 1) +
    geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "red") +
    facet_wrap(~chromosome, scales = "free_x", ncol = 5) +
    scale_color_manual(values = rep(color_palette, length.out = length(unique(eqtl_data$chromosome)))) +
    theme_minimal() +
    labs(
      title = "Manhattan Plot - cis-eQTL Analysis",
      x = "Genomic Position",
      y = "-log10 P-value",
      color = "Chromosome"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 10, face = "bold")
    )
  
  save_plot(manhattan_plot, filename, width = 15, height = 10)
  return(manhattan_plot)
}

# Function to create Venn diagram
create_venn_diagram <- function(gene_lists, filename = "venn_diagram.png") {
  # Create Venn diagram
  venn_plot <- venn.diagram(
    x = gene_lists,
    category.names = names(gene_lists),
    filename = file.path(figures_dir, filename),
    output = TRUE,
    imagetype = "png",
    height = 480,
    width = 480,
    resolution = 300,
    compression = "lzw",
    lwd = 2,
    lty = 'blank',
    fill = color_palette[1:length(gene_lists)],
    cex = 1.5,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    rotation = 1
  )
  
  cat("Created Venn diagram:", filename, "\n")
  return(venn_plot)
}

# =============================================================================
# Main Visualization Pipeline
# =============================================================================

cat("Starting visualization pipeline...\n")

# Check if data files exist and create sample visualizations
# For demonstration purposes, we'll create some example plots

# 1. Create sample correlation heatmap (if data exists)
cat("Creating sample correlation heatmap...\n")
# This would use actual data when available
# create_correlation_heatmap(sample_data, "sample_correlation_heatmap.png")

# 2. Create PCA plot (if data exists)
cat("Creating PCA plot...\n")
# This would use actual data when available
# create_pca_plot(sample_data, sample_groups, "pca_plot.png")

# 3. Create volcano plot (if DEG data exists)
cat("Creating volcano plot...\n")
# This would use actual data when available
# create_volcano_plot(deg_data, "volcano_plot.png")

# 4. Create Manhattan plot (if eQTL data exists)
cat("Creating Manhattan plot...\n")
# This would use actual data when available
# create_manhattan_plot(eqtl_data, "manhattan_plot.png")

# 5. Create Venn diagram (if multiple gene lists exist)
cat("Creating Venn diagram...\n")
# This would use actual data when available
# create_venn_diagram(gene_lists, "venn_diagram.png")

# =============================================================================
# Create Example Plots for Demonstration
# =============================================================================

cat("Creating example plots for demonstration...\n")

# Example 1: Sample correlation heatmap with simulated data
set.seed(42)
n_samples <- 20
n_genes <- 100
sample_data <- matrix(rnorm(n_samples * n_genes), nrow = n_genes, ncol = n_samples)
colnames(sample_data) <- paste0("Sample_", 1:n_samples)
rownames(sample_data) <- paste0("Gene_", 1:n_genes)

create_correlation_heatmap(sample_data, "example_correlation_heatmap.png")

# Example 2: PCA plot with simulated data
sample_groups <- rep(c("Control", "Infected"), each = n_samples/2)
create_pca_plot(sample_data, sample_groups, "example_pca_plot.png")

# Example 3: Volcano plot with simulated DEG data
deg_data <- data.frame(
  gene = paste0("Gene_", 1:1000),
  log2FoldChange = rnorm(1000, mean = 0, sd = 2),
  padj = runif(1000, 0, 1)
)
create_volcano_plot(deg_data, "example_volcano_plot.png")

# Example 4: Manhattan plot with simulated eQTL data
eqtl_data <- data.frame(
  chromosome = rep(1:22, each = 50),
  position = rep(1:50, 22) * 1000000,
  pvalue = runif(1100, 0, 1)
)
create_manhattan_plot(eqtl_data, "example_manhattan_plot.png")

# Example 5: Venn diagram with simulated gene lists
gene_lists <- list(
  "Healthy_vs_Infected" = paste0("Gene_", 1:500),
  "Young_vs_Old" = paste0("Gene_", 300:800),
  "Caucasian_vs_AfricanAmerican" = paste0("Gene_", 200:700)
)
create_venn_diagram(gene_lists, "example_venn_diagram.png")

# =============================================================================
# Create Summary Report
# =============================================================================

cat("Creating visualization summary report...\n")

# Generate summary of created plots
plots_created <- list.files(figures_dir, pattern = "*.png")
summary_file <- file.path(figures_dir, "visualization_summary.txt")

cat("=== Visualization Summary ===\n", file = summary_file)
cat("Total plots created:", length(plots_created), "\n", file = summary_file, append = TRUE)
cat("Plots directory:", figures_dir, "\n", file = summary_file, append = TRUE)
cat("\nPlots created:\n", file = summary_file, append = TRUE)

for (plot in plots_created) {
  cat("- ", plot, "\n", file = summary_file, append = TRUE)
}

cat("\nVisualization pipeline complete!\n")
cat("Created", length(plots_created), "plots in:", figures_dir, "\n")
cat("Check the visualization summary:", file.path(figures_dir, "visualization_summary.txt"), "\n")
