# Next Generation Sequencing and Bioinformatics Project

**Course:** LIFE622C - Next Generation Sequencing and Bioinformatics  
**Institution:** POSTECH (Pohang University of Science and Technology)  
**Duration:** August 2024 - December 2024  
**Faculty:** Dr. Jong Kyoung Kim, Laboratory of Computational Biology  
**Student:** Saggurthi Dinesh (BE21B032)

## Project Overview

This repository contains the reproduction and analysis of RNA-Seq data from the research paper "Host response to influenza infections in human blood: association of influenza severity with host genetics and transcriptomic response" (Schughart et al., 2024). The project focuses on identifying key molecular pathways and biomarkers for influenza infection through comprehensive bioinformatics analysis, including differential expression analysis, cis-eQTL analysis, and population-specific investigations.

## Research Objectives

1. **RNA-Seq Analysis Reproduction**: Reproduce the RNA-Seq analysis pipeline using modern bioinformatics tools
2. **Differential Expression Analysis**: Identify differentially expressed genes (DEGs) between healthy and infected samples
3. **cis-eQTL Analysis**: Investigate ethnicity-associated genetic variations
4. **Pathway Analysis**: Discover key molecular pathways involved in influenza response
5. **Biomarker Identification**: Identify potential biomarkers for influenza infection

## Tools and Technologies Used

### Core Bioinformatics Tools
- **STAR**: RNA-seq read alignment
- **Trim Galore**: Quality trimming and adapter removal
- **FastQC**: Quality control assessment
- **MultiQC**: Quality control report aggregation
- **DESeq2**: Differential expression analysis
- **Limma**: Linear models for microarray and RNA-seq data
- **Rsubread**: RNA-seq read counting

### Analysis and Visualization
- **R**: Statistical computing and graphics
- **R Markdown**: Reproducible research documents
- **ggplot2**: Data visualization
- **Venn diagrams**: Gene set comparisons
- **PCA plots**: Principal component analysis
- **Manhattan plots**: Genome-wide association results
- **Volcano plots**: Differential expression results

## Project Structure

```
ngs-bioinformatics-project/
├── data/                    # Data files
│   ├── raw/                # Raw sequencing data
│   ├── processed/          # Processed data files
│   └── reference/          # Reference genomes and annotations
├── scripts/                # Analysis scripts
│   ├── preprocessing/      # Data preprocessing scripts
│   ├── analysis/          # Analysis scripts
│   └── visualization/     # Plotting and visualization scripts
├── notebooks/              # Jupyter/R Markdown notebooks
├── results/                # Analysis results
│   ├── figures/           # Generated plots and figures
│   ├── tables/            # Statistical tables
│   └── reports/           # Analysis reports
├── workflows/              # Workflow definitions
└── docs/                   # Documentation
```

## Key Findings

### Differential Expression Analysis
- Identified significant DEGs between healthy controls and influenza-infected patients
- Discovered age-related differences in immune response (ICU young vs. ICU non-young)
- Found ethnicity-associated variations in gene expression patterns

### cis-eQTL Analysis
- Performed cis-eQTL analysis within 1Mb of transcription start sites
- Identified genetic variants associated with gene expression changes
- Explored population-specific genetic variations (Caucasian vs. African American)

### Pathway Analysis
- Identified key molecular pathways involved in influenza response
- Discovered interferon-related gene expression patterns
- Found seasonal variations in immune response

## Data Sources

The project uses data from the following sources:
- **ST1**: Combined transcriptome data (2018-2022)
- **ST2**: Coriell samples description
- **ST3**: Genotype data with ethnicity information
- **ST4**: Seasonal collection data
- **ST5**: Collection site information
- **ST6-ST8**: Healthy control vs. infected comparisons
- **ST9-ST11**: ICU patient comparisons
- **ST12**: Age-related analysis (old vs. young)
- **ST13**: cis-eQTL analysis results
- **ST14**: Ethnicity comparison results
- **ST20**: Normalized LIMMA batch-corrected data
- **ST21**: Combined numeric genotypes
- **ST26**: Genetic group 4 genotypes
- **ST27**: Transcriptome data
- **ST28**: Infected African American vs. Caucasian comparison
- **ST29-ST30**: Interferon response analysis

## Getting Started

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)
- Required R packages: DESeq2, Limma, Rsubread, ggplot2, dplyr, etc.

### Installation
1. Clone this repository
2. Install required R packages
3. Download and place data files in appropriate directories
4. Run analysis scripts in order

### Usage
1. Start with preprocessing scripts in `scripts/preprocessing/`
2. Run analysis scripts in `scripts/analysis/`
3. Generate visualizations using `scripts/visualization/`
4. View results in `results/` directory

## Results and Visualizations

### Key Figures
- **DEG Venn Diagram**: Overlapping differentially expressed genes across comparisons
- **PCA Plots**: Principal component analysis of sample relationships
- **Volcano Plots**: Differential expression significance vs. fold change
- **Manhattan Plots**: Genome-wide association results
- **Heatmaps**: Gene expression patterns across samples

### Statistical Results
- Differential expression statistics
- Pathway enrichment analysis
- cis-eQTL association results
- Population-specific genetic variations

## References

**Primary Research Paper:**
Schughart, K., Smith, A. M., Tsalik, E. L., Threlkeld, S. C., Sellers, S., Fischer, W. A. II, Schreiber, J., Lücke, E., Cornberg, M., Debarry, J., Woods, C. W., McClain, M. T., & Heise, M. (2024). Host response to influenza infections in human blood: association of influenza severity with host genetics and transcriptomic response. *Frontiers in Immunology*, 15, 1385362. https://doi.org/10.3389/fimmu.2024.1385362

---

*Last updated: December 2024*
