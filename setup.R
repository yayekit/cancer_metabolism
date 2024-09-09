### 0.1 load packages ###
# List of packages to load
packages <- c(
  "BiocManager", "readxl", "unix", "heatmaply", "pals", "scales", 
  "hrbrthemes", "viridis", "forcats", "ggrepel", "ggpubr", "grid", 
  "gridExtra", "splitstackshape", "RCurl", "magrittr", "tidyverse", 
  "data.table", "fuzzyjoin"
)

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Set memory limit (OS-specific)
if (.Platform$OS.type == "windows") {
  memory.size(max = TRUE)
} else {
  # Uncomment the following line for macOS/Linux:
  # rlimit_as(1e14, 1e14)
}

# Set global options
theme_set(theme_bw(base_size = 10))
options(scipen = 999, digits = 10)

# NB: RcppEigen is commented out to avoid potential conflicts
# library(RcppEigen)


### 0.2 links ###
# Define data source URLs with comments
data_urls <- list(
  ccle_metabolites = "https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv",  # Clear metabolomics data by the CCLE.
  ccle_names = "https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt",  # Cancer cell lines names for dataframes. Important for merging.
  ctrp = "ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip",  # CTRP dataset for dose-response data.
  gdsc_to_ccle_names = "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4E.xlsx",  # Names for overlapping with CCLE-dataframe.
  gdsc_to_ctrp_names = "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4I.xlsx"  # Names for overlapping with CTRP-dataframe.
)

# Assign URLs to individual variables
list2env(data_urls, envir = .GlobalEnv)
