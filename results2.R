results2 <- list.files("results2", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  format_results()

results2.summary <- results2 %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
  print(n = Inf) 

barplot_causal_all <- function(results) {
  
  results %>%
    barplot_causal() +
    geom_hline(aes(yintercept = auc_max), linetype = 3, color = "blue",
               data = data.frame(par.h2 = c(0.5, 0.8), auc_max = c(0.84, 0.94))) +
    facet_grid(par.dist ~ par.h2)
}

results2.summary %>%
  filter(par.model == "fancy", method %in% c("logit-simple", "PRS-max")) %>%
  barplot_causal_all() +
  labs(y = "Mean of 20 AUCs")
  
ggsave("figures/supp-AUC-logit-fancy.pdf", scale = 1/90, width = 1242, height = 951)

results2.summary %>%
  filter(par.model == "fancy", grepl("PRS", method)) %>%
  barplot_causal_all() +
  labs(y = "Mean of 20 AUCs")

ggsave("figures/supp-AUC-PRS-fancy.pdf", scale = 1/90, width = 1242, height = 951)

## Logistic regressions

compare_logit(bind_rows(results1, results2), with_method = "logit-triple")

ggsave("figures/supp-triple.pdf", scale = 1/90, width = 1400, height = 1060)

results2.summary %>%
  filter(grepl("logit", method)) %>%
  barplot_causal_one(h2 = 0.5) +
  facet_grid(par.dist ~ par.model) +
  labs(y = "Mean of 20 AUCs")

ggsave("figures/supp-AUC-triple.pdf", scale = 1/90, width = 1242, height = 951)
