# Data Directory

This directory contains the data files used in the Next Generation Sequencing and Bioinformatics project.

## Data Sources

The project uses data from the FImm (Finnish Immune Response to Influenza) study, which investigates host response to influenza infections in human blood.

## Data Files Structure

### Core Data Files

| File | Description | Size | Format |
|------|-------------|------|--------|
| ST1_transcriptome_sbst1_target_comb_SIG_2018_2019_2020_2022_181122_2a.txt | Combined transcriptome data (2018-2022) | Large | Tab-separated |
| ST2_coriell_samples_description_kls_310823.xlsx | Coriell samples description | Medium | Excel |
| ST3_genotype_target_combined_unique_SIG_geno_assethn_310823a.txt | Genotype data with ethnicity information | Large | Tab-separated |
| ST4_season_110823.xlsx | Seasonal collection data | Small | Excel |
| ST5_Collection_site_040923.xlsx | Collection site information | Small | Excel |

### Analysis Results

| File | Description | Analysis Type |
|------|-------------|---------------|
| ST6_sbst1_limma_healthy_control_vs_infected_181122.txt | Healthy control vs. infected comparison | Differential expression |
| ST7_sbst1_limma_ICUn_vsCNRTL_181122.txt | ICU non-young vs. control comparison | Differential expression |
| ST8_sbst1_limma_ICUy_vsCNRTL_181122.txt | ICU young vs. control comparison | Differential expression |
| ST9_sbst1_limma_ICUy_vs_ICUn_181122.txt | ICU young vs. ICU non-young comparison | Differential expression |
| ST10_sbst1_limma_ICUn_fVSm_231122.txt | ICU non-young vs. ICU young comparison | Differential expression |
| ST11_sbst1_limma_ICUy_fVSm_231122.txt | ICU young vs. ICU non-young comparison | Differential expression |
| ST12_sbst1_limma_old_young_060923.txt | Age-related analysis (old vs. young) | Differential expression |
| ST13_cis_eQTL_DEGs_dist1Mb_result_011222_150923.txt | cis-eQTL analysis results | Genetic association |
| ST14_DEGs_CaucVSAfrAM_ciseQTL_121023.txt | Ethnicity comparison results | Population genetics |
| ST20_sbst1_norm_LIMMA_btch_corr_SIG_comb_2022_2020_2019_2018_181122_1.txt | Normalized LIMMA batch-corrected data | Normalized expression |
| ST21_comb_numeric_genotypes_SIG_geno_ptntID_011222.7z | Combined numeric genotypes | Genotype data |
| ST26_gg4_geno_011222.7z | Genetic group 4 genotypes | Genotype data |
| ST27_tr7_transcriptome_011222.txt | Transcriptome data | Expression data |
| ST28_sbst1_limma_infAfrAm_vs_infCaucs_181122a.txt | Infected African American vs. Caucasian comparison | Population-specific |
| ST29_ifn_old_071123.txt | Interferon response in old samples | Pathway analysis |
| ST30_ifn_young_071123.txt | Interferon response in young samples | Pathway analysis |

## Data Organization

### Raw Data
- **Sequencing files**: FASTQ files from RNA-Seq experiments
- **Reference genomes**: Human genome reference files
- **Annotation files**: Gene annotation and metadata

### Processed Data
- **Expression matrices**: Normalized gene expression data
- **Quality metrics**: QC reports and statistics
- **Alignment files**: BAM/SAM files from STAR alignment

### Analysis Results
- **Differential expression**: DESeq2 and Limma results
- **Genetic associations**: eQTL analysis results
- **Pathway analysis**: Gene set enrichment results

## Data Formats

### Text Files (.txt)
- Tab-separated values
- UTF-8 encoding
- Header row with column names

### Excel Files (.xlsx)
- Multiple sheets for different data types
- Metadata and sample information
- Collection and processing details

### Compressed Files (.7z, .zip)
- Large genotype and sequence data
- Use appropriate tools for extraction
- Check file integrity after extraction

## Data Quality

### Quality Control Metrics
- **FastQC reports**: Raw sequencing quality
- **MultiQC summaries**: Aggregated QC metrics
- **Sample correlation**: Inter-sample relationships
- **PCA analysis**: Sample clustering and outliers

### Data Validation
- **Expression ranges**: Check for reasonable values
- **Missing data**: Identify and handle missing values
- **Batch effects**: Assess and correct for technical variation
- **Sample metadata**: Verify sample information consistency

## Usage Guidelines

### Data Access
1. **Raw data**: Access through secure channels only
2. **Processed data**: Available for analysis scripts
3. **Results**: Generated and stored in results directory

### Data Processing
1. **Preprocessing**: Use scripts in `scripts/preprocessing/`
2. **Analysis**: Use scripts in `scripts/analysis/`
3. **Visualization**: Use scripts in `scripts/visualization/`

### Data Storage
1. **Large files**: Store in appropriate locations
2. **Backup**: Regular backup of important data
3. **Version control**: Track data processing steps

## References

1. **Primary Research Paper**: Schughart, K., et al. (2024). Host response to influenza infections in human blood: association of influenza severity with host genetics and transcriptomic response. *Frontiers in Immunology*, 15, 1385362. https://doi.org/10.3389/fimmu.2024.1385362
2. FImm Study Protocol
3. RNA-Seq Analysis Pipeline Documentation
4. Data Processing Standards
5. Quality Control Guidelines

