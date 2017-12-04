library(tidyverse)

format_results <- function(results) {
  results %>%
    as_tibble() %>%
    mutate(
      par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                          levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
      AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2]))
    )
}

results3 <- list.files(paste0("results", 2:3), full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  filter(method != "logit-triple", par.model != "fancy") %>%
  format_results()


results3.summary <- results3 %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
  print(n = Inf)

## Make a function of this

barplot_causal <- function(results, h2) {
  
  auc_max <- `if`(h2 == 0.8, 0.94, 0.84)
  
  results %>%
    filter(par.h2 == h2) %>%
    myggplot(aes(par.causal, AUC_mean, fill = method, color = method)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_hline(yintercept = auc_max, linetype = 3, color = "blue") +
    geom_bar(stat = "identity", alpha = 0.3, position = position_dodge()) +
    geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
                  position=position_dodge(width=0.9), color = "black", width = 0.2) +
    facet_grid(par.dist ~ .) + 
    scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
                       oob = scales::rescale_none) +
    theme(strip.text.x = element_text(size = rel(2)),
          strip.text.y = element_text(size = rel(2)))+
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
         fill = "Method", color = "Method")
}

results3.summary %>%
  filter(method %in% c("logit-simple", "PRS-max")) %>%
  barplot_causal(h2 = 0.8) +
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2)

ggsave("figures/main-AUC-logit.pdf", scale = 1/90, width = 844, height = 872)

results3.summary %>%
  filter(grepl("PRS", method)) %>%
  barplot_causal(h2 = 0.8)

ggsave("figures/main-AUC-PRS.pdf", scale = 1/90, width = 844, height = 872)

results3.summary %>%
  filter(method %in% c("logit-simple", "PRS-max")) %>%
  barplot_causal(h2 = 0.5) +
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2)

ggsave("figures/supp-AUC-logit.pdf", scale = 1/90, width = 844, height = 872)

results3.summary %>%
  filter(grepl("PRS", method)) %>%
  barplot_causal(h2 = 0.5)

ggsave("figures/supp-AUC-PRS.pdf", scale = 1/90, width = 844, height = 872)


## Corr measures
results3_corr <- results3 %>%
  filter(method %in% c("logit-simple", "PRS-max")) %>%
  mutate(
    percCases10 = map_dbl(eval, ~mean(.x[order(.x[, 1], decreasing = TRUE)[1:110], 2])),
    percCases20 = map_dbl(eval, ~mean(.x[order(.x[, 1], decreasing = TRUE)[1:220], 2])),
    Parameters  = interaction(paste("h2 =", par.h2), par.dist, method, sep = " | ")
  )
  
p1 <- results3_corr %>%
  myggplot(aes(AUC, percCases10, color = Parameters)) + 
  geom_point(size = 0.6, alpha = 0.5) +
  geom_smooth(aes(linetype = Parameters), method = "lm", se = FALSE, size = 2) +
  theme(legend.key.width = unit(3.5, "line"), 
        legend.direction = "horizontal") + 
  labs(y = "Proportion of cases among top 10% scores")
p2 <- results3_corr %>%
  myggplot(aes(AUC, percCases20, color = Parameters)) + 
  geom_point(size = 0.6, alpha = 0.5) +
  geom_smooth(aes(linetype = Parameters), method = "lm", se = FALSE, size = 2) + 
  labs(y = "Proportion of cases among top 20% scores") +
  theme(legend.position = "none")

cowplot::plot_grid(p1 + theme(legend.position = "none"), p2 , 
                   ncol = 2, align = "hv", scale = 0.95) %>%
  cowplot::plot_grid(cowplot::get_legend(p1), ., rel_heights = c(0.1, 1), ncol = 1)


ggsave("figures/supp-AUC-corr.pdf", scale = 1/90, width = 1558, height = 737)

with(results3_corr, cor(AUC, percCases10))
with(results3_corr, cor(AUC, percCases20))
