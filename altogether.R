### 0.1 libraries ###
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
library(RcppEigen)
memory.size(max = TRUE) # Windows-speific memory limit
# rlimit_as(1e14, 1e14) # macOS/Linux specific
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
# "https://ctd2-data.nci.nih.gov/Public/Broad/CTRPv2.0_2015_ctd2_ExpandedDataset/CTRPv2.0_2015_ctd2_ExpandedDataset.zip"
#
gdsc_to_ccle_names_link <-
  "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4E.xlsx"
#
gdsc_to_ctrp_names_link <-
  "https://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources//Data/suppData/TableS4I.xlsx"















##### 1.1 ctrp_list #####
# Acquiring CTRP data from the web.
tempfile_ctrp <-
  tempfile() # Allocating a temporary file to write the downloaded data to.
#
download.file(
  ctrp_link,
  tempfile_ctrp,
  mode = "wb"
) # Downloading the data and writing it to the allocated file.
#
tempdir_ctrp <-
  tempdir() # Creating a temporary directory to save downloaded files to.
#
unzip(
  tempfile_ctrp,
  exdir = "tempdir_ctrp"
) # Uncompressing zipped downloaded data.
# Attaching individual uncompressed files to a list.
ctrp_names_list <-
  list.files(
    path = "tempdir_ctrp",
    pattern = "*.txt",
    full.names = TRUE,
    recursive = FALSE
  )
# Looping over the elements of the list, reading each of the elements as a TAB-delimited dataframe
ctrp_list <-
  lapply(ctrp_names_list, function(x) {
    read_delim(
      x,
      "\t",
      escape_double = FALSE,
      trim_ws = TRUE,
      col_names = TRUE
    )
  }) %>%
  purrr::compact() # Using compact() to discard elements of the list of dataframes hat are NULL or that have length zero.
#
dir.create(path = "Data/tempdir_ctrp") # Creating a permanent directory to save processed files to.
#
ctrp_list <-
  setNames(
    object = ctrp_list,
    ctrp_names_list
  ) # Provide names to the elements of the list.
#
list2env(
  ctrp_list,
  globalenv()
) # From a named list ctrp_list, assign all list components as objects into a pre-existing environment.
#
obj <-
  ls(ctrp_list) # save the elemens of the list to a new object.
#
for (i in 1:length(obj)) {
  # save every element of the list as a separate .RData file. These files could be accessed offline later from the "Data" directory.
  save(
    list = (obj[i]),
    file = paste(
      "Data/",
      obj[i],
      ".RData",
      sep = ""
    )
  )
}
# Removing every temporary and/or intermediate file and folder used for data downloading and unpacking.
rm(
  tempfile_ctrp,
  tempdir_ctrp,
  ctrp_names_list,
  ctrp_list,
  obj
)















##### 1.2 ccle_metabolomics #####
# Downloading core metabolomics data for this project.
#
ccle_names <-
  fread(ccle_names_link) %>%
  select(
    CCLE_ID,
    Name,
    tcga_code
  )















ccle_metabolomics <-
  fread(ccle_metabolites_link) %>%
  stringdist_join(
    ccle_names,
    by = "CCLE_ID",
    max_dist = 1,
    method = "lcs",
    mode = "left",
    ignore_case = TRUE
  ) %>%
  select(-CCLE_ID.y) %>%
  rename(CCLE_ID = CCLE_ID.x) %>%
  # distinct(
  #   DepMap_ID,
  #   .keep_all = TRUE
  # ) %>%
  select(
    CCLE_ID,
    Name,
    DepMap_ID,
    tcga_code,
    everything()
  ) %>%
  rename(
    ccl_name = Name
  ) %>%
  drop_na(ccl_name)
#
save(
  ccle_metabolomics,
  file = "Data/ccle_metabolomics.RData"
)
#
rm(
  ccle_names,
  ccle_metabolites_link,
  ccle_names_link
)










##### 1.3 ctrp_for_cor_long #####
ctrp_for_cor_long <-
  `tempdir_ctrp/v20.data.curves_post_qc.txt` %>%
  full_join(`tempdir_ctrp/v20.meta.per_experiment.txt`) %>%
  full_join(`tempdir_ctrp/v20.meta.per_compound.txt`) %>%
  full_join(`tempdir_ctrp/v20.meta.per_cell_line.txt`) %>%
  select(
    apparent_ec50_umol,
    area_under_curve,
    top_test_conc_umol,
    cpd_name,
    master_cpd_id,
    broad_cpd_id,
    ccl_name,
    master_ccl_id,
    growth_mode,
    # Pathway.Name,
    ccle_primary_site,
    ccle_primary_hist,
    ccle_hist_subtype_1,
    gene_symbol_of_protein_target
    # tcga_code
  ) %>%
  drop_na(ccl_name) %>%
  stringdist_join(
    ccle_metabolomics,
    by = "ccl_name",
    max_dist = 1,
    method = "lcs",
    mode = "right",
    ignore_case = TRUE
  ) %>%
  rename(ccl_name = as.factor("ccl_name.y")) %>%
  select(-c(
    # ccl_name.y,
    ccl_name.x,
    CCLE_ID,
    DepMap_ID
  )) %>%
  pivot_longer(
    !c(
      apparent_ec50_umol,
      area_under_curve,
      top_test_conc_umol,
      cpd_name,
      master_cpd_id,
      broad_cpd_id,
      ccl_name,
      master_ccl_id,
      growth_mode,
      # Pathway.Name,
      ccle_primary_site,
      ccle_primary_hist,
      ccle_hist_subtype_1,
      gene_symbol_of_protein_target,
      tcga_code
    ),
    names_to = "metabolite_name",
    names_repair = "check_unique",
    values_to = "metabolite_conc",
    values_drop_na = FALSE
  ) %>%
  group_by(
    tcga_code,
    master_cpd_id,
    metabolite_name
  ) %>%
  mutate(
    group_id = cur_group_id() # Gives a unique numeric identifier for each group
  ) %>%
  ungroup()
#
save(
  ctrp_for_cor_long,
  file = "Data/ctrp_for_cor_long.RData"
)
#
rm(
  `tempdir_ctrp/v20.data.curves_post_qc.txt`,
  `tempdir_ctrp/v20.meta.per_experiment.txt`,
  `tempdir_ctrp/v20.meta.per_compound.txt`,
  # `tempdir_ctrp/v20.meta.per_cell_line.txt`,
  `tempdir_ctrp/v20._COLUMNS.txt`,
  `tempdir_ctrp/v20._README.txt`,
  `tempdir_ctrp/v20.data.per_cpd_avg.txt`,
  `tempdir_ctrp/v20.data.per_cpd_post_qc.txt`,
  `tempdir_ctrp/v20.data.per_cpd_pre_qc.txt`,
  `tempdir_ctrp/v20.data.per_cpd_well.txt`,
  `tempdir_ctrp/v20.meta.media_comp.txt`,
  `tempdir_ctrp/v20.meta.per_assay_plate.txt`,
  ccle_metabolomics,
  ctrp_link,
  i
)










##### 2 ctrp_cor #####
ctrp_for_cor <-
  ctrp_for_cor_long %>%
  drop_na() %>%
  mutate(tcga_code = as.factor(tcga_code)) %>%
  group_by(tcga_code) %>%
  filter(n_distinct(master_ccl_id) >= 15) %>%
  ungroup() %>%
  filter(!str_detect(cpd_name, ":")) %>%
  group_by(group_id) %>%
  filter(n() >= 15) %>%
  ungroup()
#
save(
  ctrp_for_cor,
  file = "Data/ctrp_for_cor.RData"
)










library(RcppEigen)
##### 2 ctrp_cor #####
ctrp_cor <-
  setDT(ctrp_for_cor)[, .(
    cor_coef = cor.test(
      metabolite_conc,
      area_under_curve,
      method = "pearson",
      alternative = "two.sided",
      na.action = "na.exclude"
    )$estimate,
    cor_p_val = cor.test(
      metabolite_conc,
      area_under_curve,
      method = "pearson",
      alternative = "two.sided",
      na.action = "na.exclude"
    )$p.value,
    lm_p_val = as.numeric(summary(fastLm(
      cbind(1, metabolite_conc), area_under_curve
    ))$coefficients[2, 4]),
    y_intercept = fastLmPure(cbind(1, metabolite_conc), area_under_curve)$coefficients[1],
    lin_gradient = fastLmPure(cbind(1, metabolite_conc), area_under_curve)$coefficients[2],
    lin_rmse = fastLmPure(cbind(1, metabolite_conc), area_under_curve)$s,
    metabolite_conc = metabolite_conc,
    area_under_curve = area_under_curve,
    apparent_ec50_umol = apparent_ec50_umol,
    top_test_conc_umol = top_test_conc_umol,
    cpd_name = cpd_name,
    master_cpd_id = master_cpd_id,
    broad_cpd_id = broad_cpd_id,
    master_ccl_id = master_ccl_id,
    ccl_name = ccl_name,
    growth_mode = growth_mode,
    metabolite_name = metabolite_name,
    # Metabolite.ID = Metabolite.ID,
    # Pathway.Name = Pathway.Name,
    # Pathway.Subject = Pathway.Subject,
    ccle_primary_site = ccle_primary_site,
    ccle_primary_hist = ccle_primary_hist,
    ccle_hist_subtype_1 = ccle_hist_subtype_1,
    gene_symbol_of_protein_target = gene_symbol_of_protein_target,
    tcga_code = tcga_code
  ), by = "group_id"] %>%
  group_by(
    tcga_code,
    master_cpd_id
  ) %>%
  mutate(
    cor_q_val = p.adjust(
      cor_p_val,
      method = "BH"
    ),
    lm_q_val = p.adjust(
      lm_p_val,
      method = "BH"
    )
  ) %>%
  ungroup()
#
detach(
  "package:RcppEigen",
  unload = TRUE
)
#
save(
  ctrp_cor,
  file = "Data/ctrp_cor.RData"
)
#
rm(
  ctrp_for_cor_long,
  ctrp_for_cor,
)










ctrp_cor_short <-
  ctrp_cor %>%
  group_by(group_id) %>%
  distinct(
    cor_coef,
    .keep_all = TRUE
  ) %>%
  ungroup()
#
save(
  ctrp_cor_short,
  file = "Data/ctrp_cor_short.RData"
)










detach("package:dplyr", unload = TRUE)
library(plyr) # Specific methods of data binding, not availabele in dplyr
library(dplyr)
#
metabol_link <-
  "https://smpdb.ca/downloads/smpdb_metabolites.csv.zip"
#
# load("C:/Users/yemelianovskyi/Google Drive/Study/ICB/thesisProject/Data/ctrp_cor.RData")
#
##### 1.1 ctrp_list #####
# Acquiring CTRP data from the web.
tempfile_metabol <-
  tempfile() # Allocating a temporary file to write the downloaded data to.
#
getBinaryURL(
  metabol_link,
  ftp.use.epsv = FALSE,
  crlf = TRUE
) %>%
  writeBin(con = tempfile_metabol) # Downloading the data and writing it to the allocated file.
#
tempdir_metabol <-
  tempdir() # Creating a temporary directory to save downloaded files to.
#
unzip(
  tempfile_metabol,
  exdir = "tempdir_metabol"
) # Uncompressing zipped downloaded data.
#
# REPEATABLE!
#
metabol_names_list <-
  list.files(
    path = "tempdir_metabol",
    pattern = "\\.csv$",
    full.names = TRUE,
    recursive = FALSE
  )
#
metabol_names_smpdb <-
  lapply(
    metabol_names_list,
    fread
  )
#
names(metabol_names_smpdb) <-
  metabol_names_list
#
metabol_names_smpdb <-
  ldply(
    metabol_names_smpdb,
    data.frame,
    .id = "SMPDB.ID"
  )
#
save(
  metabol_names_smpdb,
  file = "Data/metabol_names_smpdb.RData"
)










metabol_names_and_pathways <-
  metabol_names_smpdb %>%
  rename(metabolite_name = Metabolite.Name) %>%
  select(
    # SMPDB.ID,
    Pathway.Name,
    Pathway.Subject,
    Metabolite.ID,
    metabolite_name
  ) %>%
  filter(
    Pathway.Subject %in% c(
      "Disease",
      # "Metabolic"
      "Drug Action",
      "Signaling",
      "Protein",
      "Physiological",
      "Drug Metabolism"
    )
  ) %>%
  group_by(
    Pathway.Subject,
    Metabolite.ID,
    metabolite_name
  ) %>%
  summarise(
    # SMPDB.ID = paste(unique(SMPDB.ID), collapse = '; '),
    Pathway.Name = paste(unique(Pathway.Name), collapse = "; ")
  ) %>%
  ungroup()
# distinct(
#   metabolite_name,
#   .keep_all = TRUE
# )
#
save(
  metabol_names_and_pathways,
  file = "Data/metabol_names_and_pathways.RData"
)










ctrp_metabolites_cor_short <-
  ctrp_cor_short %>%
  # mutate(
  #   metabolite_name = as.factor(tolower(Metabolite.Name)),
  #   .keep = "unused"
  # ) %>%
  stringdist_join(
    metabol_names_and_pathways,
    by = "metabolite_name",
    max_dist = 2,
    method = "lcs",
    mode = "left",
    ignore_case = TRUE
  ) %>%
  rename(metabolite_name = metabolite_name.x) %>%
  mutate(metabolite_name.y = NULL)
#
save(
  ctrp_metabolites_cor_short,
  file = "Data/ctrp_metabolites_cor_short.RData"
)










##### 11 ctrp_metabolites_cor_signif #####
ctrp_metabolites_cor_signif <-
  ctrp_cor %>%
  ungroup() %>%
  mutate(correlation = case_when(
    cor_coef <= -0.4 &
      cor_q_val < 0.1
    ~ "strong_negative",
    cor_coef >= 0.4 &
      cor_q_val < 0.1
    ~ "strong_positive",
    TRUE
    ~ "insignificant"
  ))
#
save(
  ctrp_metabolites_cor_signif,
  file = "Data/ctrp_metabolites_cor_signif.RData"
)










##### 11 ctrp_cor_short_signif #####
ctrp_cor_short_signif <-
  ctrp_cor_short %>%
  ungroup() %>%
  mutate(
    correlation = case_when(
      (
        cor_coef <= -0.4 &
          cor_q_val < 0.1
      )
      ~ "strong_negative",
      (
        cor_coef >= 0.4 &
          cor_q_val < 0.1
      )
      ~ "strong_positive",
      TRUE
      ~ "insignificant"
    )
  )
# filter(
#   point_colour != "grey"
# )
#
save(
  ctrp_cor_short_signif,
  file = "Data/ctrp_cor_short_signif.RData"
)










##### 11 ctrp_metabolites_cor_signif_PDF #####
ctrp_metabolites_cor_signif_PDF <-
  ctrp_metabolites_cor_signif %>%
  group_by(group_id) %>%
  mutate(
    conc_min = min(as.vector(metabolite_conc)),
    conc_max = max(as.vector(metabolite_conc)),
    conc_range = conc_max - conc_min,
    auc_mean_value = mean(as.vector(area_under_curve))
  ) %>%
  # ungroup() %>%
  arrange(
    auc_mean_value,
    desc(conc_range),
    cor_q_val,
    .by_group = TRUE
  ) %>%
  split(forcats::fct_inorder(factor(.$group_id)))
#
save(
  ctrp_metabolites_cor_signif_PDF,
  file = "Data/ctrp_metabolites_cor_signif_PDF.RData"
)







































hyper_target_prep <-
  ctrp_metabolites_cor_short %>%
  select(
    "group_id",
    "cor_coef",
    "cor_q_val",
    "master_cpd_id",
    "master_ccl_id",
    "ccl_name",
    "gene_symbol_of_protein_target",
    "tcga_code",
    "Metabolite.ID",
    "Pathway.Name"
  ) %>%
  drop_na() %>%
  # Separating all the targets from the "gene_symbol_of_protein_target" column, if there is more than one target per cell
  # distinct(
  #   group_id,
  #   .keep_all = TRUE
  # ) %>%
  separate_rows(
    gene_symbol_of_protein_target,
    sep = ";"
  ) %>%
  # mutate(
  #   gene_symbol_of_protein_target = as.factor(gene_symbol_of_protein_target),
  #   Pathway.Name = as.factor(Pathway.Name)
  # ) %>%
  distinct(
    gene_symbol_of_protein_target,
    master_ccl_id,
    master_cpd_id,
    .keep_all = TRUE
  ) %>%
  group_by(tcga_code) %>%
  # How many times is this target aimed by *this particular* drug, producing *an SCCLP response*?
  add_tally(name = "ALL_ALL_TARGETS") %>%
  # mutate(ALL_ALL_TARGETS = n_distinct(master_cpd_id)) %>%
  ungroup() %>%
  # How many times is this target aimed by *this particular* drug, producing *an SCCLP response*?
  group_by(
    tcga_code,
    gene_symbol_of_protein_target
  ) %>%
  add_tally(name = "ALL_PER_TARGET") %>%
  # mutate(ALL_PER_TARGET = n_distinct(master_cpd_id)) %>%
  ungroup()
#
save(
  hyper_target_prep,
  file = "Data/hyper_target_prep.RData"
)















hyper_target_positive <-
  hyper_target_prep %>%
  # How many times is this target aimed by *this particular* drug, producing *any kind of response*?
  # How many times is this target aimed by *any* drug, producing *any kind of response*?
  group_by(
    tcga_code,
    cor_q_val < 0.1,
    cor_coef > 0
  ) %>%
  add_tally(name = "SIG_ALL_TARGETS") %>%
  # mutate(SIG_ALL_TARGETS = n_distinct(master_cpd_id)) %>%
  ungroup() %>%
  dplyr::filter(
    cor_q_val < 0.1,
    cor_coef > 0
  ) %>%
  # distinct(
  #   gene_symbol_of_protein_target,
  #   master_ccl_id,
  #   master_cpd_id,
  #   .keep_all = TRUE
  # ) %>%
  group_by(
    tcga_code,
    gene_symbol_of_protein_target
  ) %>%
  add_tally(name = "SIG_PER_TARGET") %>%
  # mutate(SIG_PER_TARGET = n_distinct(master_cpd_id)) %>%
  ungroup() %>%
  mutate(
    # Here the hypergeometric score is being calculated. Explaining it with the "marble sack" example:
    P_TARGET = phyper(
      # q = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls
      q = SIG_PER_TARGET - 1,
      # m = the number of white balls in the urn.
      m = SIG_ALL_TARGETS,
      # n = the number of black balls in the urn
      n = ALL_ALL_TARGETS - SIG_ALL_TARGETS,
      # k = the number of balls drawn from the urn
      k = ALL_PER_TARGET,
      # Setting lower.tail as FALSE, because if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
      lower.tail = FALSE,
      log.p = FALSE
    )
  ) %>%
  filter(SIG_PER_TARGET > 3) %>%
  group_by(tcga_code) %>%
  # Adjusting the hypergeometric p-value using the BH method (https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
  mutate(
    ADJ_P_TARGET = p.adjust(
      p = P_TARGET,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  # To make sure that we do not include duplicated rows to the results
  distinct(
    ADJ_P_TARGET,
    gene_symbol_of_protein_target,
    .keep_all = TRUE
  ) %>%
  # Leaving only significantly enriched (P < 0.05) drugs
  dplyr::filter(ADJ_P_TARGET < 0.1)
#
save(
  hyper_target_positive,
  file = "Data/hyper_target_positive.RData"
)















hyper_target_negative <-
  hyper_target_prep %>%
  # How many times is this target aimed by *this particular* drug, producing *any kind of response*?
  # How many times is this target aimed by *any* drug, producing *any kind of response*?
  group_by(
    tcga_code,
    cor_q_val < 0.1,
    cor_coef < 0
  ) %>%
  add_tally(name = "SIG_ALL_TARGETS") %>%
  # mutate(SIG_ALL_TARGETS = n_distinct(master_cpd_id)) %>%
  ungroup() %>%
  dplyr::filter(
    cor_q_val < 0.1,
    cor_coef < 0
  ) %>%
  # distinct(
  #   gene_symbol_of_protein_target,
  #   master_ccl_id,
  #   master_cpd_id,
  #   .keep_all = TRUE
  # ) %>%
  group_by(
    tcga_code,
    gene_symbol_of_protein_target
  ) %>%
  add_tally(name = "SIG_PER_TARGET") %>%
  # mutate(SIG_PER_TARGET = n_distinct(master_cpd_id)) %>%
  ungroup() %>%
  mutate(
    # Here the hypergeometric score is being calculated. Explaining it with the "marble sack" example:
    P_TARGET = phyper(
      # q = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls
      q = SIG_PER_TARGET - 1,
      # m = the number of white balls in the urn.
      m = SIG_ALL_TARGETS,
      # n = the number of black balls in the urn
      n = ALL_ALL_TARGETS - SIG_ALL_TARGETS,
      # k = the number of balls drawn from the urn
      k = ALL_PER_TARGET,
      # Setting lower.tail as FALSE, because if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
      lower.tail = FALSE,
      log.p = FALSE
    )
  ) %>%
  filter(SIG_PER_TARGET > 3) %>%
  group_by(tcga_code) %>%
  # Adjusting the hypergeometric p-value using the BH method (https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
  mutate(
    ADJ_P_TARGET = p.adjust(
      p = P_TARGET,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  # To make sure that we do not include duplicated rows to the results
  distinct(
    ADJ_P_TARGET,
    gene_symbol_of_protein_target,
    .keep_all = TRUE
  ) %>%
  # Leaving only significantly enriched (P < 0.05) drugs
  dplyr::filter(ADJ_P_TARGET < 0.1)
#
save(
  hyper_target_negative,
  file = "Data/hyper_target_negative.RData"
)















hyper_target <-
  hyper_target_positive %>%
  bind_rows(hyper_target_negative) %>%
  mutate(
    sign = case_when(
      cor_coef < 0 ~ "-",
      cor_coef > 0 ~ "+"
    )
  ) %>%
  # Ordering the sresults nicely
  select(
    tcga_code,
    gene_symbol_of_protein_target,
    # Pathway.Name,
    sign,
    cor_coef,
    cor_q_val,
    P_TARGET,
    ADJ_P_TARGET,
    SIG_PER_TARGET,
    ALL_PER_TARGET,
    SIG_ALL_TARGETS,
    ALL_ALL_TARGETS,
    ccl_name,
    group_id
    # everything()
  ) %>%
  # Arranging the results so the lowest adjusted hypergeometric p-values are on the top, and the highest are on the bottom of the data frame
  arrange(
    ADJ_P_TARGET,
    gene_symbol_of_protein_target
  )
#
save(
  hyper_target,
  file = "Data/hyper_target.RData"
)


































### 7.1 hyper_pathway_prep ###
library(tidyverse)
library(stringr)
library(ggrepel)
memory.limit(size = 56000)
theme_set(theme_bw(base_size = 18))
options(scipen = 0, digits = 5)
#
hyper_pathway_prep <-
  ctrp_metabolites_cor_short %>%
  select(
    "group_id",
    "cor_coef",
    "cor_q_val",
    "master_cpd_id",
    "master_ccl_id",
    "ccl_name",
    "gene_symbol_of_protein_target",
    "tcga_code",
    "Metabolite.ID",
    "Pathway.Name"
  ) %>%
  drop_na() %>%
  # Separating all the targets from the "gene_symbol_of_protein_target" column, if there is more than one target per cell
  # distinct(
  #   group_id,
  #   .keep_all = TRUE
  # ) %>%
  separate_rows(
    gene_symbol_of_protein_target,
    sep = ";"
  ) %>%
  # mutate(
  #   gene_symbol_of_protein_target = as.factor(gene_symbol_of_protein_target),
  #   Pathway.Name = as.factor(Pathway.Name)
  # ) %>%
  distinct(
    tcga_code,
    master_cpd_id,
    Pathway.Name,
    Metabolite.ID,
    .keep_all = TRUE
  ) %>%
  group_by(tcga_code) %>%
  # How many times is this target aimed by *this particular* drug, producing *an SCCLP response*?
  # add_tally(name = "ALL_ALL_TARGETS") %>%
  mutate(ALL_ALL_TARGETS = n_distinct(Metabolite.ID)) %>%
  ungroup() %>%
  # How many times is this target aimed by *this particular* drug, producing *an SCCLP response*?
  group_by(
    tcga_code,
    Pathway.Name
  ) %>%
  # add_tally(name = "ALL_PER_TARGET") %>%
  mutate(ALL_PER_TARGET = n_distinct(Metabolite.ID)) %>%
  ungroup()
#
save(
  hyper_pathway_prep,
  file = "Data/hyper_pathway_prep.RData"
)















hyper_pathway_positive <-
  hyper_pathway_prep %>%
  # How many times is this target aimed by *this particular* drug, producing *any kind of response*?
  # How many times is this target aimed by *any* drug, producing *any kind of response*?
  group_by(
    tcga_code,
    cor_q_val < 0.1,
    cor_coef > 0
  ) %>%
  # add_tally(name = "SIG_ALL_TARGETS") %>%
  mutate(SIG_ALL_TARGETS = n_distinct(Metabolite.ID)) %>%
  ungroup() %>%
  dplyr::filter(
    cor_q_val < 0.1,
    cor_coef > 0
  ) %>%
  # distinct(
  #   Pathway.Name,
  #   master_ccl_id,
  #   master_cpd_id,
  #   .keep_all = TRUE
  # ) %>%
  group_by(
    tcga_code,
    Pathway.Name
  ) %>%
  # add_tally(name = "SIG_PER_TARGET") %>%
  mutate(SIG_PER_TARGET = n_distinct(Metabolite.ID)) %>%
  ungroup() %>%
  mutate(
    # Here the hypergeometric score is being calculated. Explaining it with the "marble sack" example:
    P_PATHWAY = phyper(
      # q = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls
      q = SIG_PER_TARGET - 1,
      # m = the number of white balls in the urn.
      m = SIG_ALL_TARGETS,
      # n = the number of black balls in the urn
      n = ALL_ALL_TARGETS - SIG_ALL_TARGETS,
      # k = the number of balls drawn from the urn
      k = ALL_PER_TARGET,
      # Setting lower.tail as FALSE, because if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
      lower.tail = FALSE,
      log.p = FALSE
    )
  ) %>%
  filter(SIG_PER_TARGET > 3) %>%
  group_by(tcga_code) %>%
  # Adjusting the hypergeometric p-value using the BH method (https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
  mutate(
    ADJ_P_PATHWAY = p.adjust(
      p = P_PATHWAY,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  # To make sure that we do not include duplicated rows to the results
  distinct(
    ADJ_P_PATHWAY,
    Pathway.Name,
    .keep_all = TRUE
  ) %>%
  # Leaving only significantly enriched (P < 0.05) drugs
  dplyr::filter(ADJ_P_PATHWAY < 0.1)
#
save(
  hyper_pathway_positive,
  file = "Data/hyper_pathway_positive.RData"
)















hyper_pathway_negative <-
  hyper_pathway_prep %>%
  # How many times is this target aimed by *this particular* drug, producing *any kind of response*?
  # How many times is this target aimed by *any* drug, producing *any kind of response*?
  group_by(
    tcga_code,
    cor_q_val < 0.1,
    cor_coef < 0
  ) %>%
  # add_tally(name = "SIG_ALL_TARGETS") %>%
  mutate(SIG_ALL_TARGETS = n_distinct(Metabolite.ID)) %>%
  ungroup() %>%
  dplyr::filter(
    cor_q_val < 0.1,
    cor_coef < 0
  ) %>%
  # distinct(
  #   Pathway.Name,
  #   master_ccl_id,
  #   master_cpd_id,
  #   .keep_all = TRUE
  # ) %>%
  group_by(
    tcga_code,
    Pathway.Name
  ) %>%
  # add_tally(name = "SIG_PER_TARGET") %>%
  mutate(SIG_PER_TARGET = n_distinct(Metabolite.ID)) %>%
  ungroup() %>%
  mutate(
    # Here the hypergeometric score is being calculated. Explaining it with the "marble sack" example:
    P_PATHWAY = phyper(
      # q = vector of quantiles representing the number of white balls drawn without replacement from an urn which contains both black and white balls
      q = SIG_PER_TARGET - 1,
      # m = the number of white balls in the urn.
      m = SIG_ALL_TARGETS,
      # n = the number of black balls in the urn
      n = ALL_ALL_TARGETS - SIG_ALL_TARGETS,
      # k = the number of balls drawn from the urn
      k = ALL_PER_TARGET,
      # Setting lower.tail as FALSE, because if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x]
      lower.tail = FALSE,
      log.p = FALSE
    )
  ) %>%
  filter(SIG_PER_TARGET > 3) %>%
  group_by(tcga_code) %>%
  # Adjusting the hypergeometric p-value using the BH method (https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html)
  mutate(
    ADJ_P_PATHWAY = p.adjust(
      p = P_PATHWAY,
      method = "BH"
    )
  ) %>%
  ungroup() %>%
  # To make sure that we do not include duplicated rows to the results
  distinct(
    ADJ_P_PATHWAY,
    Pathway.Name,
    .keep_all = TRUE
  ) %>%
  # Leaving only significantly enriched (P < 0.05) drugs
  dplyr::filter(ADJ_P_PATHWAY < 0.1)
#
save(
  hyper_pathway_negative,
  file = "Data/hyper_pathway_negative.RData"
)















hyper_pathway <-
  hyper_pathway_positive %>%
  bind_rows(hyper_pathway_negative) %>%
  mutate(
    sign = case_when(
      cor_coef < 0 ~ "-",
      cor_coef > 0 ~ "+"
    )
  ) %>%
  # Ordering the sresults nicely
  select(
    tcga_code,
    Pathway.Name,
    sign,
    cor_coef,
    cor_q_val,
    P_PATHWAY,
    ADJ_P_PATHWAY,
    SIG_PER_TARGET,
    ALL_PER_TARGET,
    SIG_ALL_TARGETS,
    ALL_ALL_TARGETS,
    ccl_name,
    group_id
    # everything()
  ) %>%
  # Arranging the results so the lowest adjusted hypergeometric p-values are on the top, and the highest are on the bottom of the data frame
  arrange(
    ADJ_P_PATHWAY,
    Pathway.Name
  )
#
save(
  hyper_pathway,
  file = "Data/hyper_pathway.RData"
)















hyper_both <-
  hyper_target %>%
  merge(
    hyper_pathway,
    by = c(
      "tcga_code",
      "sign"
      # "ccl_name"
      # "Pathway.Name"
    )
  ) %>%
  select(
    tcga_code,
    gene_symbol_of_protein_target,
    Pathway.Name,
    sign,
    ADJ_P_TARGET,
    ADJ_P_PATHWAY,
    P_TARGET,
    P_PATHWAY,
    cor_coef.x,
    cor_coef.y,
    cor_q_val.x,
    cor_q_val.y
  ) %>%
  distinct(across(everything(), .keep_all = TRUE)) %>%
  group_by(across(c(
    tcga_code,
    gene_symbol_of_protein_target,
    sign,
    ADJ_P_TARGET,
    ADJ_P_PATHWAY,
    cor_coef.x,
    cor_coef.y,
    cor_q_val.x,
    cor_q_val.y
  ))) %>%
  summarise(
    # gene_symbol_of_protein_target = paste(unique(gene_symbol_of_protein_target), collapse = ";",
    Pathway.Name = paste(unique(Pathway.Name), collapse = "; ")
    # Pathway.Name = paste(unique(Pathway.Name), collapse = "; ")
  ) %>%
  ungroup()
#
save(
  hyper_both,
  file = "Data/hyper_both.RData"
)





























volcano_by_tcga_code <-
  hyper_both %>%
  mutate(text = fct_reorder(tcga_code, cor_coef.x)) %>%
  ggplot(aes(
    x = cor_coef.x,
    y = -log(cor_q_val.x, 10),
    color = text,
    fill = text
  ),
  xlim = c(-1, 1)
  # ylim = c(0, 1)
  ) +
  geom_point(alpha = 1, size = 2) +
  # geom_vline(
  #   aes(xintercept = 0.5),
  #   colour = "black",
  #   linetype = "twodash",
  #   size = 0.5,
  #   show.legend = FALSE
  # ) +
  geom_vline(
    aes(xintercept = 0),
    colour = "black",
    linetype = "solid",
    size = 0.5,
    show.legend = FALSE
  ) +
  # geom_vline(
  #   aes(xintercept = -0.5),
  #   colour = "black",
  #   linetype = "twodash",
  #   size = 0.5,
  #   show.legend = FALSE
  # ) +
  geom_hline(
    aes(yintercept = -log10(0.1)),
    colour = "black",
    linetype = "twodash",
    size = 0.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = as.vector(polychrome(36))) +
  scale_fill_manual(values = as.vector(polychrome(36))) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 11, face = "bold"),
    legend.title = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(colour = "black", size = 12)
  ) +
  # xlim = c(0, max(-log(ctrp_cor$cor_q_val, 10))) +
  xlab("Pearson correlation score (on a scale from -1 to +1") +
  ylab("Negative logarithm base 10 from FDR corrected Pearson correlation p-value (q-value)") +
  facet_grid(~text) +
  theme_bw()
#
save(
  volcano_by_tcga_code,
  file = "Data/volcano_by_tcga_code.RData"
)
#
volcano_by_tcga_code










# hyper_both <-
#   hyper_both %>%
#   group_by(ADJ_P_PATHWAY) %>%
#   mutate(gene_symbol_of_protein_target = paste0(gene_symbol_of_protein_target, collapse = ";")) %>%
#   mutate_at(vars(Pathway.Name), funs(replace(Pathway.Name, duplicated(Pathway.Name), NA)))
# #
# save(
#   hyper_both,
#   file = "Data/hyper_both.RData"
# )










##### 11 summary_logq #####
summary_logq <-
  ggplot(
    hyper_both,
    aes(
      x = -log(ADJ_P_TARGET, 10),
      y = -log(ADJ_P_PATHWAY, 10)
    ),
    # xlim = c(-1, 1),
    ylim = c(-20, 70)
  ) +
  geom_point(aes(colour = tcga_code),
    alpha = 0.75,
    size = 4
  ) +
  geom_vline(
    aes(xintercept = 0.1),
    colour = "black",
    linetype = "twodash",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = -0.1),
    colour = "black",
    linetype = "twodash",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_hline(
    aes(yintercept = 0.1),
    colour = "black",
    linetype = "twodash",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_hline(
    aes(yintercept = -log10(ADJ_P_PATHWAY)),
    linetype = "dotted",
    size = 0.25
  ) +
  geom_label_repel(
    aes(-log10(ADJ_P_TARGET),
      -log10(ADJ_P_PATHWAY),
      label = Pathway.Name
    ),
    size = 2.5,
    fontface = "bold",
    nudge_x = 0.15,
    nudge_y = 0.5,
    segment.size = 0.5,
    segment.curvature = -0.1,
    segment.alpha = 1,
    segment.color = "black",
    force = 1,
    force_pull = 0,
    max.iter = 10000000,
    max.overlaps = 25,
    min.segment.length = 0,
    arrow = arrow(length = unit(0.01, "npc")),
    direction = "both"
  ) +
  geom_vline(aes(xintercept = -log(ADJ_P_TARGET, 10)), linetype = "dotted", size = 0.25) +
  geom_text(aes(
    -log10(ADJ_P_TARGET),
    -15,
    label = gene_symbol_of_protein_target
  ),
  size = 2.25,
  vjust = -0.5,
  angle = 90
  ) +
  # geom_label_repel(
  #   aes(
  #     x = -log(SIG_PER_TARGET.x, 10),
  #     y = -log(SIG_PER_TARGET.y, 10),
  #     label = paste0(
  #       "gene_symbol_of_protein_target = ",
  #       gene_symbol_of_protein_target,
  #       "\n",
  #       "Pathway.Name = ",
  #       Pathway.Name
  #     )
  #   ),
  #   fontface = "bold",
  #   segment.size = 1,
  #   segment.alpha = 0.8,
  #   segment.color = "grey",
  #   force = 1,
  #   force_pull = 0,
  #   max.iter = 100000,
  #   max.overlaps = 25,
  #   direction = "both",
  #   hjust = 0.5,
  #   vjust = 1.5,
  #   size = 2.25
  # ) +
  scale_colour_manual(values = as.vector(jet(7))) +
  xlab("Negative log10 of the Hypergeometric Adjusted Score of the Protein Target") +
  ylab("Negative log10 of the Hypergeometric Adjusted Score of the Metabolic Pathway")
#
# save(
#   summary_logq,
#   file = "Data/summary_logq.RData"
# )
#
summary_logq



















hist(tstats, freq = FALSE, col = "gray", breaks = 20)







options(scipen = 0, digits = 5)
#
test_metabolites_cor_short_hypergeom_combined_Target_PDF <-
  test_metabolites_cor_CCLE_NAMES %>%
  semi_join(
    (hyper_both %>%
      select(group_id.x) %>%
      rename(group_id = group_id.x)),
    by = "group_id"
  ) %>%
  group_by(group_id) %>%
  mutate(
    conc_min = min(as.vector(metabolite_conc)),
    conc_max = max(as.vector(metabolite_conc)),
    conc_range = conc_max - conc_min,
    auc_mean_value = mean(as.vector(area_under_curve))
  ) %>%
  ungroup() %>%
  arrange(
    auc_mean_value,
    desc(conc_range),
    cor_q_val
    # .by_group = TRUE
  ) %>%
  distinct(.keep_all = TRUE) %>%
  split(forcats::fct_inorder(factor(.$group_id)))
#
save(
  test_metabolites_cor_short_hypergeom_combined_Target_PDF,
  file = "Data/test_metabolites_cor_short_hypergeom_combined_Target_PDF.RData"
)























options(scipen = 0, digits = 5)
##### 11 p_ctrp_signif_dose_response_by_ccle_primary_site.pdf #####
# Flush all the results to one PDF-file
pdf(
  "Plots/p_ctrp_signif_dose_response_by_ccle_primary_site.pdf",
  onefile = TRUE,
  height = 20,
  width = 20
)
#
for (i in 1:length(test_metabolites_cor_short_hypergeom_combined_Target_PDF))
{
  df <- as.data.frame(test_metabolites_cor_short_hypergeom_combined_Target_PDF[[i]])
  #
  cell_line_broad <-
    paste("BROAD CCL ID: ", df$ccl_name)
  #
  primary_site_ccle <-
    paste("Primary site: ", df$tcga_code)
  #
  compound_id_broad <-
    paste("BROAD compound ID: ", df$broad_cpd_id)
  #
  compound_name_broad <-
    paste("Compound name: ", df$cpd_name)
  #
  metab_name <-
    paste("Metabolite name: ", df$metabolite_name)
  #
  area_under_curve_value <-
    paste("App. EC50: ", df$area_under_curve)
  #
  y_intercept_value <-
    paste("Y-intercept: ", df$y_intercept)
  #
  lin_gradient_value <-
    paste("Linear regression gradient:", df$lin_gradient)
  #
  lin_rmse_value <-
    paste("Linear RMSE: ", df$lin_rmse)
  #
  lm_q_val_value <-
    paste("Linear fit p-value: ", df$lm_q_val)
  #
  cor_q_val_value <-
    paste("Corr. p. value: ", df$cor_q_val)
  #
  cor_coef_value <-
    paste("Corr. coeff. value: ", df$cor_coef)
  #
  gene_symbol_of_protein_target_value <-
    paste("Gene of the protein target: ", df$gene_symbol_of_protein_target)
  #
  chart <-
    ggplot(df) +
    theme(
      axis.text.x = element_text(
        face = "bold",
        color = "#993333",
        size = 21
      ),
      axis.text.y = element_text(
        face = "bold",
        color = "#993333",
        size = 21,
        angle = 45
      ),
      axis.title.x = element_text(size = 21, face = "bold"),
      axis.title.y = element_text(size = 21, face = "bold"),
      plot.title = element_text(
        color = "black",
        size = 21,
        face = "bold"
      ),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)
    ) +
    theme(legend.position = "right") +
    guides(color = guide_colourbar(direction = "vertical")) +
    xlab("Metabolite concentration, µM/ml") +
    ylab("Area Under the Curve (AUC)") +
    geom_point(
      aes(
        x = metabolite_conc,
        y = area_under_curve,
        colour = cor_coef
      ),
      size = 9,
      alpha = 1
    ) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "cor_coef",
      limits = c(-1, 1)
    ) +
    geom_abline(
      intercept = df$y_intercept,
      slope = df$lin_gradient,
      colour = "green",
      size = 1.5,
      linetype = "dashed"
    ) +
    geom_label_repel(
      aes(
        x = metabolite_conc,
        y = area_under_curve,
        label = paste0(
          # "AUC = ",
          # area_under_curve,
          # "\n",
          # "CCL name = ",
          ccl_name
        )
      ),
      fontface = "bold",
      segment.size = 0.25,
      segment.alpha = 0.8,
      segment.color = "grey",
      force = 1,
      force_pull = 0,
      max.iter = 100000,
      max.overlaps = 100,
      direction = "both",
      # hjust = 0.5,
      # vjust = -1.5,
      size = 5
    ) +
    ggtitle(
      paste0(
        primary_site_ccle,
        "\n",
        metab_name,
        "\n",
        compound_id_broad,
        "\n",
        compound_name_broad,
        "\n",
        cor_coef_value,
        "\n",
        cor_q_val_value,
        "\n",
        lin_rmse_value,
        "\n",
        lm_q_val_value,
        "\n",
        gene_symbol_of_protein_target_value
      )
    )
  #
  print(chart)
}
#
dev.off()































































# Basic barplot
barplot_tcga <-
  ctrp_for_cor %>%
  group_by(tcga_code) %>%
  distinct(
    master_ccl_id,
    .keep_all = TRUE
  ) %>%
  mutate(n_distinct_ccl = n_distinct(master_ccl_id)) %>%
  ungroup() %>%
  distinct(
    tcga_code,
    .keep_all = TRUE
  ) %>%
  ggplot(., aes(x = n_distinct_ccl, y = tcga_code)) +
  geom_bar(stat = "identity")
#
barplot_tcga




















##### 3 histogram_by_tcga_code #####
histogram_by_tcga_code <-
  ctrp_cor_short_signif %>%
  # filter(area_under_curve < 100000) %>%
  mutate(text = fct_reorder(
    tcga_code,
    area_under_curve
  )) %>%
  ggplot(aes(
    x = area_under_curve,
    color = tcga_code,
    fill = tcga_code
  )) +
  geom_histogram(
    alpha = 0.6,
    bins = 100,
    binwidth = 0.1
  ) +
  scale_colour_manual(values = as.vector(polychrome(35))) +
  scale_fill_manual(values = as.vector(polychrome(35))) +
  theme(
    legend.position = "none",
    panel.spacing = unit(
      0.1,
      "lines"
    ),
    strip.text.x = element_text(
      size = 14,
      face = "bold"
    )
  ) +
  xlab("") +
  ylab("Number of responses with such AUC value, #") +
  facet_wrap(~text) +
  theme_bw()
#
histogram_by_tcga_code




















##### 11 volcano_by_tcga_code #####
volcano_by_tcga_code <-
  ctrp_cor_short_signif %>%
  mutate(
    tcga_cancer_code = fct_reorder(
      tcga_code,
      cor_coef
    )
  ) %>%
  ggplot(aes(
    x = cor_coef,
    y = -log10(cor_q_val),
    color = correlation,
    fill = correlation
  ),
  xlim = c(-1, 1),
  # ylim = c(0, 1)
  ) +
  geom_point(alpha = 1, size = 1.5) +
  geom_vline(
    aes(xintercept = 0.5),
    colour = "black",
    linetype = "solid",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = 0),
    colour = "black",
    linetype = "dotted",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = -0.5),
    colour = "black",
    linetype = "solid",
    size = 0.5,
    show.legend = FALSE
  ) +
  geom_hline(
    aes(yintercept = -log10(0.1)),
    colour = "black",
    linetype = "solid",
    size = 0.5,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "strong_negative" = "blue",
    "strong_positive" = "red",
    "insignificant" = "grey"
  )) +
  scale_fill_manual(values = c(
    "strong_negative" = "blue",
    "strong_positive" = "red",
    "insignificant" = "grey"
  )) +
  theme_ipsum() +
  theme(
    legend.position = "none",
    panel.spacing = unit(
      0.1,
      "lines"
    ),
    strip.text.x = element_text(
      size = 11,
      face = "bold"
    ),
    legend.title = element_text(
      colour = "black",
      size = 12,
      face = "bold"
    ),
    legend.text = element_text(
      colour = "black",
      size = 12
    )
  ) +
  # xlim = c(0, max(-log(ctrp_cor$cor_q_val, 10))) +
  xlab("Pearson correlation score (on a scale from -1 to +1)") +
  ylab("Negative logarithm base 10 from FDR corrected Pearson correlation p-value (q-value)") +
  facet_wrap(~tcga_cancer_code) +
  theme_bw()
#
volcano_by_tcga_code




















##### 11 volcano_summary #####
volcano_summary <- ggplot(
  ctrp_cor_short_signif,
  aes(
    x = cor_coef,
    y = -log10(cor_q_val),
  ),
) +
  geom_point(aes(
    colour = correlation,
    fill = tcga_code
  ),
  shape = 21,
  alpha = 0.75,
  stroke = 2,
  size = 5
  ) +
  geom_vline(
    aes(xintercept = 0.4),
    colour = "black",
    linetype = "twodash",
    size = 1,
    show.legend = FALSE
  ) +
  geom_vline(
    aes(xintercept = -0.4),
    colour = "black",
    linetype = "twodash",
    size = 1,
    show.legend = FALSE
  ) +
  geom_hline(
    aes(yintercept = -log10(0.1)),
    colour = "black",
    linetype = "twodash",
    size = 1,
    show.legend = FALSE
  ) +
  scale_colour_manual(values = c(
    "strong_negative" = "blue",
    "strong_positive" = "red",
    "insignificant" = "grey"
  )) +
  scale_fill_manual(values = as.vector(polychrome(35))) +
  theme_bw()
#
volcano_summary




















##### 0.1 libraries #####
library(tidyverse)
library(data.table)
library(readxl)
library(RCurl) # Obtaining the CCLE metabolic data from the web.
library(heatmaply)
library(pals)
library(scales)
library(hrbrthemes)
library(viridis)
library(forcats)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(grid)
library(gridExtra)
library(splitstackshape)
theme_set(theme_bw(base_size = 10))
options(scipen = 0, digits = 10)
memory.limit(size = 56000)




hypergeom_group_id <-
  hyper_both %>%
  select(
    tcga_code,
    gene_symbol_of_protein_target,
    Pathway.Name
  )






##### 11 test_metabolites_base #####
test_metabolites_base <-
  ctrp_metabolites_cor_short %>%
  separate_rows(
    gene_symbol_of_protein_target,
    sep = ";",
    convert = TRUE
  ) %>%
  # separate_rows(
  #   Pathway.Name,
  #   sep = "; ",
  #   convert = TRUE
  # ) %>%
  as.data.frame() %>%
  right_join(hypergeom_group_id) %>%
  select(
    tcga_code,
    metabolite_name,
    Pathway.Name,
    gene_symbol_of_protein_target,
    cpd_name,
    cor_coef
  ) %>%
  split(.$tcga_code) %>%
  Filter(nrow, .) %>%
  # Remove empty elements of the list
  imap(
    ~ pivot_wider(
      .x,
      id_cols = c(
        Pathway.Name,
        metabolite_name
      ),
      names_from = c(
        cpd_name,
        gene_symbol_of_protein_target
      ),
      values_from = cor_coef
    ) %>%
      mutate(
        across(
          where(is.character),
          as.factor
        )
      ) %>%
      as.data.frame() %>%
      column_to_rownames(var = "metabolite_name")
  )








test_metabolites_heatmap_COLS_ANNOTATION <-
  test_metabolites_base %>%
  imap(
    ~ select(
      .x,
      !Pathway.Name
    ) %>%
      t() %>%
      as.data.frame() %>%
      rownames_to_column(var = "drug_gene") %>%
      select(drug_gene) %>%
      separate(
        col = "drug_gene",
        into = c(
          "drug",
          "gene"
        ),
        sep = "_",
        remove = TRUE,
        convert = TRUE,
        fill = "left"
      ) %>%
      pivot_wider(
        id_cols = drug,
        names_from = gene,
        values_from = gene
      ) %>%
      column_to_rownames(var = "drug") %>%
      mutate_all(
        funs(
          case_when(
            is.na(.) ~ FALSE,
            TRUE ~ TRUE
          )
        )
      )
  )










test_metabolites_heatmap_ROWS_ANNOTATION <-
  test_metabolites_base %>%
  imap(
    ~ select(
      .x,
      "Pathway.Name"
    ) %>%
      rename(is = "Pathway.Name") %>%
      cSplit_e(
        split.col = "is",
        type = "character",
        drop = TRUE,
        fill = 0
      ) %>%
      as.data.frame() %>%
      mutate(
        across(
          is.numeric,
          as.logical
        )
      )
  )










test_metabolites_heatmap_prep <-
  test_metabolites_base %>%
  imap(
    ~ arrange(
      .x,
      Pathway.Name
    ) %>%
      select(
        # .x,
        !Pathway.Name
      ) %>%
      .[!duplicated(as.list(.))] %>%
      rename_all(
        ~ str_replace_all(
          string = .,
          pattern = "\\_[^_]+$",
          replacement = ""
        )
      )
  )










test_metabolites_heatmap_plot <-
  test_metabolites_heatmap_prep %>%
  imap(
    ~ heatmaply_cor(
      .x,
      # colors = ocean.balance(1000),
      na.value = "darkslategrey",
      dist_method = "maximum",
      hclust_method = "ward.D2",
      dendrogram = "both",
      k_col = NA,
      k_row = NA,
      # Rowv = FALSE,
      row_dend_left = FALSE,
      hide_colorbar = FALSE,
      # colorbar_xanchor = "middle",
      # colorbar_yanchor = "top",
      col_side_colors = test_metabolites_heatmap_COLS_ANNOTATION[[.y]],
      row_side_colors = test_metabolites_heatmap_ROWS_ANNOTATION[[.y]],
      na.rm = TRUE,
      xlab = paste0(
        "cpd_name",
        "\n",
        .y
      ),
      ylab = "metabolite_name",
      key.title = "cor_coef",
      seriate = "mean",
      plot_method = "ggplot",
      fontsize_row = 10,
      fontsize_col = 10,
      side_color_colorbar_len = 1
    )
  )
#
test_metabolites_heatmap_plot

























##### 11 p_ctrp_signif_dose_response_by_ccle_primary_site.pdf #####
# Flush all the results to one PDF-file
pdf(
  "Plots/p_ctrp_signif_dose_response_by_ccle_primary_site.pdf",
  onefile = TRUE,
  height = 20,
  width = 20
)
#
for (i in 1:length(ctrp_metabolites_cor_signif_PDF))
{
  df <- as.data.frame(ctrp_metabolites_cor_signif_PDF[[i]])
  #
  cell_line_broad <-
    paste("BROAD CCL ID: ", df$master_ccl_id)
  #
  primary_site_ccle <-
    paste("Primary site: ", df$tcga_code)
  #
  compound_id_broad <-
    paste("BROAD compound ID: ", df$broad_cpd_id)
  #
  compound_name_broad <-
    paste("Compound name: ", df$cpd_name)
  #
  metab_name <-
    paste("Metabolite name: ", df$metabolite_name)
  #
  area_under_curve_value <-
    paste("App. EC50: ", df$area_under_curve)
  #
  y_intercept_value <-
    paste("Y-intercept: ", df$y_intercept)
  #
  lin_gradient_value <-
    paste("Linear regression gradient:", df$lin_gradient)
  #
  lin_rmse_value <-
    paste("Linear RMSE: ", df$lin_rmse)
  #
  lm_q_val_value <-
    paste("Linear fit p-value: ", df$lm_q_val)
  #
  cor_q_val_value <-
    paste("Corr. p. value: ", df$cor_q_val)
  #
  cor_coef_value <-
    paste("Corr. coeff. value: ", df$cor_coef)
  #
  chart <-
    ggplot(df) +
    theme(
      axis.text.x = element_text(
        face = "bold",
        color = "#993333",
        size = 21
      ),
      axis.text.y = element_text(
        face = "bold",
        color = "#993333",
        size = 21,
        angle = 45
      ),
      axis.title.x = element_text(size = 21, face = "bold"),
      axis.title.y = element_text(size = 21, face = "bold"),
      plot.title = element_text(
        color = "black",
        size = 21,
        face = "bold"
      ),
      legend.title = element_text(size = 14, face = "bold"),
      legend.text = element_text(size = 14)
    ) +
    theme(legend.position = "right") +
    guides(color = guide_colourbar(direction = "vertical")) +
    xlab("Metabolite concentration, µM/ml") +
    ylab("Area Under the Curve (AUC)") +
    geom_point(
      aes(
        x = metabolite_conc,
        y = area_under_curve,
        colour = cor_coef
      ),
      size = 9,
      alpha = 1
    ) +
    scale_color_gradient2(
      low = "blue",
      mid = "white",
      high = "red",
      midpoint = 0,
      name = "cor_coef",
      limits = c(-1, 1)
    ) +
    geom_abline(
      intercept = df$y_intercept,
      slope = df$lin_gradient,
      colour = "green",
      size = 1.5,
      linetype = "dashed"
    ) +
    geom_label_repel(
      aes(
        x = metabolite_conc,
        y = area_under_curve,
        label = paste0(
          # "AUC = ",
          # area_under_curve,
          # "\n",
          "CCL ID = ",
          master_ccl_id
        )
      ),
      fontface = "bold",
      segment.size = 0.25,
      segment.alpha = 0.8,
      segment.color = "grey",
      force = 1,
      force_pull = 0,
      max.iter = 100000,
      max.overlaps = 1,
      direction = "both",
      hjust = 0.5,
      vjust = -1.5,
      size = 4.5
    ) +
    ggtitle(
      paste0(
        primary_site_ccle,
        "\n",
        metab_name,
        "\n",
        compound_id_broad,
        "\n",
        compound_name_broad,
        "\n",
        cor_coef_value,
        "\n",
        cor_q_val_value,
        "\n",
        lin_rmse_value,
        "\n",
        lm_q_val_value
      )
    )
  #
  print(chart)
}
#
dev.off()
