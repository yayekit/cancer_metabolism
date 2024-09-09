source("setup.R")
load("Data/ctrp_for_cor_long.RData")

### 1.3 ctrp_for_cor_long ###
ctrp_for_cor_long <- `tempdir_ctrp/v20.data.curves_post_qc.txt` %>%
  full_join(`tempdir_ctrp/v20.meta.per_experiment.txt`) %>%
  full_join(`tempdir_ctrp/v20.meta.per_compound.txt`) %>%
  full_join(`tempdir_ctrp/v20.meta.per_cell_line.txt`) %>%
  select(master_cpd_id, cpd_name, broad_cpd_id, apparent_ec50_umol, 
         auc, ccle_primary_site, ccle_name, master_ccl_id, 
         experiment_id, master_cell_id) %>%
  left_join(ccle_metabolomics) %>%
  filter(!is.na(metabolite_name)) %>%
  mutate(group_id = paste(master_cpd_id, metabolite_name, sep = "_"))

# Save processed data
save(ctrp_for_cor_long, file = "Data/ctrp_for_cor_long.RData")

# Clean up
rm(list = setdiff(ls(), "ctrp_for_cor_long"))