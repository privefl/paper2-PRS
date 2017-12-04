library(tidyverse)

results1 <- list.files("results1", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  format_results()

results_ttrees <- results1 %>%
  filter(method %in% c("logit-simple", "T-Trees")) %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(
    AUC_mean = mean(AUC), 
    AUC_boot = boot(AUC, 1e5, mean), 
    N = n(),
    nb_preds_mean = mean(nb.preds),
    timing_mean = mean(timing)
  ) %>%
  print(n = Inf)

p1 <- barplot_causal(results_ttrees, h2 = 0.8) +
  facet_grid(par.model ~ par.dist) +
  labs(y = "Mean of 5 AUCs")

plot_ttrees <- function(results, y, ylab = y) {
  
  results %>%
    filter(method %in% c("logit-simple", "T-Trees")) %>%
    myggplot() +
    geom_boxplot(aes_string("par.causal", y, color = "method", 
                            fill = "method"), alpha = 0.3) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(par.model ~ par.dist) +
    theme(strip.text.x = element_text(size = rel(2)),
          strip.text.y = element_text(size = rel(2))) +
    labs(x = "Causal SNPs (number and location)", y = ylab, 
         fill = "Method", color = "Method") +
    theme(legend.position = "none")
}


p2 <- cowplot::plot_grid(
  plot_ttrees(results1, "nb.preds", "Number of predictors"),
  plot_ttrees(results1, "timing", "Execution time (in seconds)"),
  ncol = 1, align = "hv", scale = 0.95, labels = LETTERS[2:3], label_size = 15
)
  
p12 <- cowplot::plot_grid(
  p1 + theme(legend.position = "none"), p2, 
  ncol = 2, scale = 0.95, rel_widths = c(6, 4), labels = c("A", "")
)

cowplot::plot_grid(
  cowplot::get_legend(p1 + theme(legend.direction = "horizontal")), p12, 
  rel_heights = c(0.1, 1), ncol = 1
)

ggsave("figures/supp-ttrees.pdf", scale = 1/90, width = 1400, height = 1060)
