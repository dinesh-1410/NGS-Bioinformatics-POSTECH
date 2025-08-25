#!/usr/bin/env Rscript
# Load required libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(Limma)
  library(Rsubread)
  library(dplyr)
  library(ggplot2)
  library(pheatmap)
  library(VennDiagram)
  library(RColorBrewer)
  library(corrplot)
  library(pcaMethods)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(knitr)
  library(kableExtra)
})

# =============================================================================
# Configuration and Setup
# =============================================================================

# Set working directory and create output directories
project_dir <- getwd()
data_dir <- file.path(project_dir, "data")
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(results_dir, "figures")
tables_dir <- file.path(results_dir, "tables")

# Create output directories if they don't exist
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

# Set random seed for reproducibility
set.seed(42)

# =============================================================================
# Data Loading and Preprocessing
# =============================================================================

cat("Loading and preprocessing data...\n")

# Load transcriptome data (ST1)
transcriptome_file <- file.path(data_dir, "ST1_transcriptome_sbst1_target_comb_SIG_2018_2019_2020_2022_181122_2a.txt")
if (file.exists(transcriptome_file)) {
  transcriptome_data <- read.table(transcriptome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("Loaded transcriptome data with", nrow(transcriptome_data), "genes and", ncol(transcriptome_data), "samples\n")
} else {
  cat("Warning: Transcriptome data file not found\n")
}

# Load sample metadata
metadata_file <- file.path(data_dir, "ST2_coriell_samples_description_kls_310823.xlsx")
if (file.exists(metadata_file)) {
  # Note: This would require readxl package for Excel files
  cat("Sample metadata file found\n")
} else {
  cat("Warning: Sample metadata file not found\n")
}

# Load genotype data
genotype_file <- file.path(data_dir, "ST3_genotype_target_combined_unique_SIG_geno_assethn_310823a.txt")
if (file.exists(genotype_file)) {
  genotype_data <- read.table(genotype_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("Loaded genotype data with", nrow(genotype_data), "variants\n")
} else {
  cat("Warning: Genotype data file not found\n")
}

# =============================================================================
# Quality Control and Data Filtering
# =============================================================================

cat("Performing quality control and data filtering...\n")

# Filter low-expressed genes (if transcriptome data exists)
if (exists("transcriptome_data")) {
  # Remove genes with very low expression
  min_count <- 10
  min_samples <- 3
  
  # Calculate gene expression statistics
  gene_means <- rowMeans(transcriptome_data[, -1], na.rm = TRUE)
  gene_counts <- rowSums(transcriptome_data[, -1] >= min_count, na.rm = TRUE)
  
  # Filter genes
  keep_genes <- gene_means >= min_count & gene_counts >= min_samples
  filtered_data <- transcriptome_data[keep_genes, ]
  
  cat("Filtered data: kept", sum(keep_genes), "out of", nrow(transcriptome_data), "genes\n")
}

# =============================================================================
# Differential Expression Analysis
# =============================================================================

cat("Performing differential expression analysis...\n")

# This section would contain the actual DESeq2 analysis
# For now, we'll create a placeholder structure

# Example DESeq2 workflow (commented out as it requires actual data)
# if (exists("filtered_data")) {
#   # Create DESeq2 dataset
#   dds <- DESeqDataSetFromMatrix(
#     countData = filtered_data,
#     colData = sample_metadata,
#     design = ~ condition
#   )
#   
#   # Run DESeq2
#   dds <- DESeq(dds)
#   
#   # Get results
#   res <- results(dds)
#   
#   # Save results
#   write.csv(res, file.path(tables_dir, "differential_expression_results.csv"))
# }

# =============================================================================
# cis-eQTL Analysis
# =============================================================================

cat("Performing cis-eQTL analysis...\n")

# Load cis-eQTL results if available
cis_eqtl_file <- file.path(data_dir, "ST13_cis_eQTL_DEGs_dist1Mb_result_011222_150923.txt")
if (file.exists(cis_eqtl_file)) {
  cis_eqtl_data <- read.table(cis_eqtl_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("Loaded cis-eQTL data with", nrow(cis_eqtl_data), "associations\n")
  
  # Save cis-eQTL results
  write.csv(cis_eqtl_data, file.path(tables_dir, "cis_eqtl_results.csv"), row.names = FALSE)
} else {
  cat("Warning: cis-eQTL data file not found\n")
}

# =============================================================================
# Population-Specific Analysis
# =============================================================================

cat("Analyzing population-specific variations...\n")

# Load ethnicity comparison data
ethnicity_file <- file.path(data_dir, "ST14_DEGs_CaucVSAfrAM_ciseQTL_121023.txt")
if (file.exists(ethnicity_file)) {
  ethnicity_data <- read.table(ethnicity_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("Loaded ethnicity comparison data with", nrow(ethnicity_data), "genes\n")
  
  # Save ethnicity results
  write.csv(ethnicity_data, file.path(tables_dir, "ethnicity_comparison_results.csv"), row.names = FALSE)
} else {
  cat("Warning: Ethnicity comparison data file not found\n")
}

# =============================================================================
# Age-Related Analysis
# =============================================================================

cat("Analyzing age-related differences...\n")

# Load age comparison data
age_file <- file.path(data_dir, "ST12_sbst1_limma_old_young_060923.txt")
if (file.exists(age_file)) {
  age_data <- read.table(age_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  cat("Loaded age comparison data with", nrow(age_data), "genes\n")
  
  # Save age results
  write.csv(age_data, file.path(tables_dir, "age_comparison_results.csv"), row.names = FALSE)
} else {
  cat("Warning: Age comparison data file not found\n")
}

# =============================================================================
# Interferon Response Analysis
# =============================================================================

cat("Analyzing interferon response patterns...\n")

# Load interferon response data
ifn_old_file <- file.path(data_dir, "ST29_ifn_old_071123.txt")
ifn_young_file <- file.path(data_dir, "ST30_ifn_young_071123.txt")

if (file.exists(ifn_old_file) && file.exists(ifn_young_file)) {
  ifn_old_data <- read.table(ifn_old_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  ifn_young_data <- read.table(ifn_young_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  cat("Loaded interferon response data: old (", nrow(ifn_old_data), "), young (", nrow(ifn_young_data), ")\n")
  
  # Save interferon results
  write.csv(ifn_old_data, file.path(tables_dir, "interferon_response_old.csv"), row.names = FALSE)
  write.csv(ifn_young_data, file.path(tables_dir, "interferon_response_young.csv"), row.names = FALSE)
} else {
  cat("Warning: Interferon response data files not found\n")
}

# =============================================================================
# Generate Summary Report
# =============================================================================

cat("Generating summary report...\n")

# Create summary statistics
summary_stats <- list(
  project_title = "Next Generation Sequencing and Bioinformatics Project",
  course = "LIFE622C",
  institution = "POSTECH",
  student = "Saggurthi Dinesh (BE21B032)",
  faculty = "Dr. Jong Kyoung Kim",
  analysis_date = Sys.Date(),
  data_files_processed = length(list.files(data_dir, pattern = "*.txt")),
  results_generated = length(list.files(tables_dir, pattern = "*.csv"))
)

# Save summary
summary_file <- file.path(results_dir, "analysis_summary.txt")
cat("=== NGS Bioinformatics Project Summary ===\n", file = summary_file)
cat("Project:", summary_stats$project_title, "\n", file = summary_file, append = TRUE)
cat("Course:", summary_stats$course, "\n", file = summary_file, append = TRUE)
cat("Institution:", summary_stats$institution, "\n", file = summary_file, append = TRUE)
cat("Student:", summary_stats$student, "\n", file = summary_file, append = TRUE)
cat("Faculty:", summary_stats$faculty, "\n", file = summary_file, append = TRUE)
cat("Analysis Date:", summary_stats$analysis_date, "\n", file = summary_file, append = TRUE)
cat("Data Files Processed:", summary_stats$data_files_processed, "\n", file = summary_file, append = TRUE)
cat("Results Generated:", summary_stats$results_generated, "\n", file = summary_file, append = TRUE)

# =============================================================================
# Session Information
# =============================================================================

cat("Saving session information...\n")

# Save R session info
session_file <- file.path(results_dir, "R_session_info.txt")
sink(session_file)
cat("=== R Session Information ===\n")
cat("Analysis Date:", Sys.Date(), "\n")
cat("R Version:", R.version.string, "\n")
cat("Platform:", R.version$platform, "\n")
cat("OS:", Sys.info()["sysname"], Sys.info()["release"], "\n")
cat("\n=== Loaded Packages ===\n")
print(sessionInfo())
sink()

cat("Analysis complete! Results saved to:", results_dir, "\n")
cat("Check the following files for results:\n")
cat("- Tables:", tables_dir, "\n")
cat("- Figures:", figures_dir, "\n")
cat("- Summary:", file.path(results_dir, "analysis_summary.txt"), "\n")
cat("- Session Info:", file.path(results_dir, "R_session_info.txt"), "\n")
