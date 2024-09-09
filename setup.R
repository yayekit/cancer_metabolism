### 0.1 load packages ###
library(BiocManager)
library(readxl)
library(unix)
library(heatmaply)
library(pals)
library(scales)
library(hrbrthemes)
library(viridis)
library(forcats)
library(ggrepel)
library(ggpubr)
library(grid)
library(gridExtra)
library(splitstackshape)
library(RCurl) # Data acquisition from the web
library(magrittr)
library(tidyverse)
library(data.table)
library(fuzzyjoin)
# library(RcppEigen) This package should probably not be attached before it's actually needed, because its functions may mask base packages functions and cause errors (see) 
memory.size(max = TRUE) # Windows-speific memory limit
# rlimit_as(1e14, 1e14) # macOS/Linux specific memory limit
theme_set(theme_bw(base_size = 10))
options(scipen = 999, digits = 10)














### 0.2 links ###
# Clear metabolomics data by the CCLE.
ccle_metabolites_link <-
  "https://data.broadinstitute.org/ccle/CCLE_metabolomics_20190502.csv"
# Cancer cell lines names for dataframes. Important for merging.
ccle_names_link <-
  "https://data.broadinstitute.org/ccle/Cell_lines_annotations_20181226.txt"
# CTRP dataset for dose-response data.
ctrp_link <-
  "ftp://caftpd.nci.nih.gov/pub/OCG-DCC/CTD2/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
# Names for overlapping with CCLE-dataframe.
gdsc_to_ccle_names_link <-
  "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4E.xlsx"
# Names for overlapping with CTRP-dataframe.
gdsc_to_ctrp_names_link <-
  "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4I.xlsx"
