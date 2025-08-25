# Next Generation Sequencing and Bioinformatics Project

## Project Overview

**Course:** LIFE622C - Next Generation Sequencing and Bioinformatics  
**Institution:** POSTECH (Pohang University of Science and Technology)  
**Duration:** August 2024 - December 2024  
**Faculty:** Dr. Jong Kyoung Kim, Laboratory of Computational Biology  
**Student:** Saggurthi Dinesh (BE21B032)

## Repository Structure

```
ngs-bioinformatics-project/
├── README.md                           # Main project documentation
├── requirements.txt                    # Required R packages and tools
├── install.R                          # R package installation script
├── PROJECT_SUMMARY.md                 # This comprehensive summary
├── .gitignore                         # Git ignore patterns
├── data/                              # Data files and documentation
│   └── README.md                      # Data directory documentation
├── scripts/                           # Analysis scripts
│   ├── analysis/                      # Main analysis scripts
│   │   └── main_analysis.R           # Primary analysis pipeline
│   ├── preprocessing/                 # Data preprocessing scripts
│   └── visualization/                 # Plotting and visualization scripts
│       └── create_plots.R            # Comprehensive visualization pipeline
├── notebooks/                         # Jupyter/R Markdown notebooks
├── workflows/                         # Workflow definitions
│   └── main_workflow.Rmd             # Complete analysis workflow
├── results/                           # Analysis results and outputs
│   ├── figures/                       # Generated plots and visualizations
│   ├── tables/                        # Statistical tables and summaries
│   ├── reports/                       # Analysis reports and documentation
│   └── README.md                      # Results directory documentation
└── docs/                              # Additional documentation
```

## Research Objectives

This project aims to reproduce and extend RNA-Seq analysis from the research paper "Host response to influenza infections in human blood: association of influenza severity with host genetics and transcriptomic response" (Schughart et al., 2024) with the following objectives:

1. **RNA-Seq Analysis Reproduction**: Reproduce the RNA-Seq analysis pipeline using modern bioinformatics tools
2. **Differential Expression Analysis**: Identify differentially expressed genes (DEGs) between healthy and infected samples
3. **cis-eQTL Analysis**: Investigate ethnicity-associated genetic variations
4. **Pathway Analysis**: Discover key molecular pathways involved in influenza response
5. **Biomarker Identification**: Identify potential biomarkers for influenza infection

## Data Sources

The project utilizes comprehensive data from the FImm (Finnish Immune Response to Influenza) study:

### Core Data Files
- **ST1**: Combined transcriptome data (2018-2022)
- **ST2**: Coriell samples description
- **ST3**: Genotype data with ethnicity information
- **ST4**: Seasonal collection data
- **ST5**: Collection site information

### Analysis Results
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

## Tools and Technologies

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

## Analysis Pipeline

### 1. Data Preprocessing
- Quality control assessment
- Data filtering and normalization
- Batch effect correction
- Sample metadata integration

### 2. Quality Control Analysis
- Sample correlation analysis
- Principal component analysis
- Outlier detection
- Quality metrics assessment

### 3. Differential Expression Analysis
- DESeq2 workflow implementation
- Statistical testing and multiple testing correction
- Effect size calculation
- Significance thresholding

### 4. cis-eQTL Analysis
- Genetic variant association testing
- Distance-based filtering
- Statistical significance assessment
- Population-specific analysis

### 5. Population and Demographics Analysis
- Ethnicity-based comparisons
- Age-related differences
- Gender-specific effects
- Clinical variable integration

### 6. Pathway and Functional Analysis
- Gene set enrichment analysis
- Pathway mapping
- Functional annotation
- Network analysis

### 7. Visualization and Reporting
- Publication-quality figures
- Interactive visualizations
- Comprehensive reporting
- Reproducible documentation

## Key Features

### Comprehensive Documentation
- Detailed README files for each directory
- Code documentation and comments
- Workflow descriptions and examples
- Installation and usage instructions

### Reproducible Research
- Version-controlled analysis scripts
- Parameter documentation
- Environment specifications
- Result tracking and validation

### Professional Standards
- Industry-standard directory structure
- Best practices for bioinformatics projects
- Quality control procedures
- Data management protocols

### Educational Value
- Clear code structure and organization
- Comprehensive examples and tutorials
- Step-by-step workflow documentation
- Best practice demonstrations

## Getting Started

### Prerequisites
- R (version 4.0 or higher)
- RStudio (recommended)
- Git for version control
- Basic knowledge of bioinformatics and R programming

### Installation
1. Clone this repository
2. Run `install.R` to install required R packages
3. Review the documentation in each directory
4. Start with the main workflow document

### Usage
1. **Data Preparation**: Place data files in appropriate directories
2. **Package Installation**: Run the installation script
3. **Analysis Execution**: Run scripts in the recommended order
4. **Result Review**: Check generated outputs and reports

## Expected Outcomes

### Analysis Results
- Differential expression statistics
- Genetic association findings
- Pathway enrichment results
- Population-specific variations

### Visualizations
- Quality control plots
- Analysis result figures
- Publication-ready graphics
- Interactive visualizations

### Documentation
- Complete analysis pipeline
- Reproducible workflows
- Method descriptions
- Result interpretations

## Quality Assurance

### Code Quality
- Consistent coding standards
- Error handling and validation
- Performance optimization
- Documentation and comments

### Data Quality
- Quality control metrics
- Data validation procedures
- Outlier detection
- Batch effect assessment

### Result Validation
- Statistical validation
- Biological plausibility
- Cross-validation approaches
- Reproducibility checks

## Future Enhancements

### Technical Improvements
- Machine learning integration
- Cloud computing support
- Database integration
- API development

### Analysis Extensions
- Multi-omics integration
- Longitudinal analysis
- Clinical validation
- Biomarker discovery

### Educational Features
- Interactive tutorials
- Video demonstrations
- Case studies
- Assessment tools

## Contributing

This is a course project repository. For questions or suggestions, please contact:
- **Course Instructor**: Dr. Jong Kyoung Kim
- **Student**: Saggurthi Dinesh (BE21B032)
- **Institution**: POSTECH

## License

This project is for educational purposes as part of the LIFE622C course at POSTECH.

## References

**Primary Research Paper:**
Schughart, K., Smith, A. M., Tsalik, E. L., Threlkeld, S. C., Sellers, S., Fischer, W. A. II, Schreiber, J., Lücke, E., Cornberg, M., Debarry, J., Woods, C. W., McClain, M. T., & Heise, M. (2024). Host response to influenza infections in human blood: association of influenza severity with host genetics and transcriptomic response. *Frontiers in Immunology*, 15, 1385362. https://doi.org/10.3389/fimmu.2024.1385362

## Acknowledgments

- **Dr. Jong Kyoung Kim** for guidance and supervision
- **POSTECH** for providing the course and resources
- **Dr. Klaus Schughart and colleagues** for the original research that inspired this reproduction

## Contact Information

**Student:** Saggurthi Dinesh  
**Email:** [Your Email]  
**Institution:** POSTECH - Pohang University of Science and Technology  
**Course:** LIFE622C - Next Generation Sequencing and Bioinformatics  
**Faculty:** Dr. Jong Kyoung Kim, Laboratory of Computational Biology

---

## Quick Start Guide

1. **Explore the Repository**: Start with `README.md` for project overview
2. **Install Dependencies**: Run `install.R` to set up R packages
3. **Review Data**: Check `data/README.md` for data structure
4. **Run Analysis**: Execute scripts in `scripts/analysis/`
5. **Generate Plots**: Use `scripts/visualization/create_plots.R`
6. **View Results**: Check the `results/` directory for outputs
7. **Follow Workflow**: Use `workflows/main_workflow.Rmd` for complete pipeline

## File Descriptions

- **`README.md`**: Main project documentation and overview
- **`requirements.txt`**: List of required R packages and tools
- **`install.R`**: Automated R package installation script
- **`main_analysis.R`**: Primary analysis pipeline script
- **`create_plots.R`**: Comprehensive visualization pipeline
- **`main_workflow.Rmd`**: Complete analysis workflow document
- **`.gitignore`**: Git ignore patterns for large data files
- **Directory READMEs**: Detailed documentation for each folder

This repository provides a complete, professional-grade bioinformatics project structure that can serve as a template for future research projects and educational purposes.