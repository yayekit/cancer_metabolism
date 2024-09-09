# Metabolomics and Drug Response Analysis in Cancer Cell Lines

This R project analyzes the relationship between metabolomics data and drug responses in cancer cell lines using various datasets, including CCLE (Cancer Cell Line Encyclopedia) and CTRP (Cancer Therapeutics Response Portal).

## Features

- Data acquisition and preprocessing from multiple sources
- Correlation analysis between metabolite concentrations and drug responses
- Hypergeometric enrichment analysis for protein targets and metabolic pathways
- Visualization of results using various plot types (heatmaps, volcano plots, etc.)

## Prerequisites

Before running this project, ensure you have R installed (version 4.0.0 or higher) along with the following R packages:

```R
install.packages(c("tidyverse", "data.table", "readxl", "RCurl", "heatmaply", "pals", "scales", "hrbrthemes", "viridis", "forcats", "EnhancedVolcano", "ComplexHeatmap", "ggrepel", "ggpubr", "grid", "gridExtra", "splitstackshape", "fuzzyjoin", "plyr"))

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("RCppEigen")
```

## Project Structure

The project consists of a single R script (`altogether.R`) that performs the following main steps:

1. Data acquisition and preprocessing
2. Correlation analysis
3. Hypergeometric enrichment analysis
4. Visualization of results

## Usage

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/metabolomics-drug-response-analysis.git
   cd metabolomics-drug-response-analysis
   ```

2. Open the `altogether.R` script in your R environment (e.g., RStudio).

3. Run the script section by section, following the comments that describe each step.

## Main Analyses

### 1. Data Acquisition and Preprocessing

- Downloads and processes CCLE metabolomics data
- Acquires and processes CTRP drug response data
- Merges datasets and performs initial data cleaning

### 2. Correlation Analysis

- Calculates correlations between metabolite concentrations and drug responses
- Filters for significant correlations

### 3. Hypergeometric Enrichment Analysis

- Performs enrichment analysis for protein targets and metabolic pathways
- Identifies significantly enriched targets and pathways

### 4. Visualization

- Generates various plots, including:
  - Volcano plots by TCGA cancer type
  - Heatmaps of correlations between metabolites and drugs
  - Scatter plots of drug responses vs. metabolite concentrations

## Output

The script generates several output files, including:

- Processed data files (saved as .RData files in the `Data/` directory)
- Visualization plots (saved as PDF files in the `Plots/` directory)

## Customization

You can customize the analysis by modifying parameters such as correlation thresholds, p-value cutoffs, and visualization settings within the script.

## Contributing

Contributions to improve the analysis pipeline are welcome. Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This project uses data from the Cancer Cell Line Encyclopedia (CCLE) and the Cancer Therapeutics Response Portal (CTRP).
- Various R packages and libraries are used for data manipulation and visualization. Please cite them appropriately in your research.
