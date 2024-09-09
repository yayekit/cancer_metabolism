# Load required libraries and data
source("setup.R")
load("Data/ctrp_for_cor_long.RData")

# Calculate correlations and summarize data
ctrp_cor <- ctrp_for_cor_long %>%
  group_by(group_id) %>%
  summarise(
    cor_pearson = cor.test(apparent_ec50_umol, metabolite_conc, method = "pearson")$estimate,
    cor_spearman = cor.test(apparent_ec50_umol, metabolite_conc, method = "spearman")$estimate,
    p_value_pearson = cor.test(apparent_ec50_umol, metabolite_conc, method = "pearson")$p.value,
    p_value_spearman = cor.test(apparent_ec50_umol, metabolite_conc, method = "spearman")$p.value,
    n = n(),
    master_cpd_id = first(master_cpd_id),
    cpd_name = first(cpd_name),
    broad_cpd_id = first(broad_cpd_id),
    metabolite_name = first(metabolite_name)
  ) %>%
  ungroup()

# Filter and adjust p-values
ctrp_cor_short <- ctrp_cor %>%
  filter(n >= 5) %>%
  mutate(
    fdr_pearson = p.adjust(p_value_pearson, method = "fdr"),
    fdr_spearman = p.adjust(p_value_spearman, method = "fdr")
  )

# Save intermediate results
save(ctrp_cor_short, file = "Data/ctrp_cor_short.RData")

# Identify significant correlations
ctrp_metabolites_cor_signif <- ctrp_cor_short %>%
  filter(fdr_pearson < 0.05 | fdr_spearman < 0.05)

# Create final dataset with significant correlations
ctrp_cor_short_signif <- ctrp_cor_short %>%
  filter(group_id %in% ctrp_metabolites_cor_signif$group_id)

# Save final results
save(ctrp_cor_short_signif, file = "Data/ctrp_cor_short_signif.RData")

# Clean up workspace
rm(list = setdiff(ls(), c("ctrp_cor", "ctrp_cor_short", "ctrp_cor_short_signif")))
