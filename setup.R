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


