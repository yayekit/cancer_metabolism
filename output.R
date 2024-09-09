source("setup.R")
load("Data/ctrp_cor_short_signif.RData")
source("visualization.R")

# Generate summary statistics
summary_stats <- ctrp_cor_short_signif %>%
  summarise(
    mean_cor_pearson = mean(cor_pearson, na.rm = TRUE),
    median_cor_pearson = median(cor_pearson, na.rm = TRUE),
    mean_cor_spearman = mean(cor_spearman, na.rm = TRUE),
    median_cor_spearman = median(cor_spearman, na.rm = TRUE),
    total_significant = n()
  )

# Write summary to file
write.csv(summary_stats, "Output/summary_statistics.csv", row.names = FALSE)

# Generate PDF report
pdf("Output/analysis_report.pdf", width = 11, height = 8.5)
print(volcano_pearson)
print(volcano_spearman)
dev.off()

# Export significant correlations
write.csv(ctrp_cor_short_signif, "Output/significant_correlations.csv", row.names = FALSE)

# Clean up
rm(list = ls())
