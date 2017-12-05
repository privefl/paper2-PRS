library(tidyverse)

results1 <- list.files("results1", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  format_results()

results_ttrees <- results1 %>%
  filter(method %in% c("logit-simple", "T-Trees")) %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
  print(n = Inf)



compare_logit <- function(results, with_method) {
  
  results <- filter(
    results, method %in% c("logit-simple", with_method), par.h2 == 0.8)
  
  results.summary <- results %>%
    group_by_at(c(vars(starts_with("par")), "method")) %>%
    summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
    print(n = Inf)
  
  nsimu <- results.summary$N[[1]]
  stopifnot(all(results.summary$N == nsimu))
  
  plot_extra <- function(y, ylab = y) {
    
    results %>%
      myggplot() +
      geom_boxplot(aes_string("par.causal", y, color = "method", 
                              fill = "method"), alpha = 0.3) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_grid(par.model ~ par.dist) +
      labs(x = "Causal SNPs (number and location)", y = ylab, 
           fill = "Method", color = "Method") +
      theme(legend.position = "none") +
      scale_fill_manual(values = methods.color) +
      scale_color_manual(values = methods.color) +
      scale_y_continuous(limits = c(0, NA))
  }
  
  p1 <- results.summary %>%
    barplot_causal_one(h2 = 0.8) +
    facet_grid(par.model ~ par.dist) +
    labs(y = sprintf("Mean of %d AUCs", nsimu))
  
  p2 <- cowplot::plot_grid(
    plot_extra("nb.preds", "Number of predictors"),
    plot_extra("timing", "Execution time (in seconds)"),
    ncol = 1, align = "hv", scale = 0.95, 
    labels = LETTERS[2:3], label_size = 15
  )
  
  p12 <- cowplot::plot_grid(
    p1 + theme(legend.position = "none"), p2, 
    ncol = 2, scale = 0.95, rel_widths = c(6, 4), 
    labels = c("A", ""), label_size = 15
  )
  
  cowplot::plot_grid(
    cowplot::get_legend(p1 + theme(legend.direction = "horizontal")), p12, 
    rel_heights = c(0.1, 1), ncol = 1
  )
}

compare_logit(results1, with_method = "T-Trees")

ggsave("figures/supp-ttrees.pdf", scale = 1/90, width = 1400, height = 1060)
