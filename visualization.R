source("setup.R")
load("Data/ctrp_cor_short_signif.RData")

# Plotting functions
plot_volcano <- function(data, x, y, color, title) {
  ggplot(data, aes(x = {{x}}, y = -log10({{y}}), color = abs({{color}}))) +
    geom_point(alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red") +
    labs(title = title, x = "Correlation", y = "-log10(p-value)") +
    theme_minimal()
}

# Create plots
volcano_pearson <- plot_volcano(ctrp_cor_short_signif, 
                                cor_pearson, 
                                p_value_pearson, 
                                cor_pearson, 
                                "Volcano Plot (Pearson)")

volcano_spearman <- plot_volcano(ctrp_cor_short_signif, 
                                 cor_spearman, 
                                 p_value_spearman, 
                                 cor_spearman, 
                                 "Volcano Plot (Spearman)")

# Save plots
ggsave("Plots/volcano_pearson.png", volcano_pearson, width = 10, height = 8)
ggsave("Plots/volcano_spearman.png", volcano_spearman, width = 10, height = 8)

# Clean up
rm(list = setdiff(ls(), c("plot_volcano", "volcano_pearson", "volcano_spearman")))
