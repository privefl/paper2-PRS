library(tidyverse)
library(bigsnpr)
library(pkg.paper.PRS)

corr <- readRDS("backingfiles/corr2.rds")
system.time(
  results2 <- list.files("results2", full.names = TRUE) %>%
    read_format_results(corr)
)  # 198 sec for 1636 of ~1MB

print(results2, width = Inf)

# saveRDS(results2, "results2_all.rds")

####

results2 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC, n = 1e4), N = n()) %>%
  print(n = Inf) -> results2.summary

barplot_causal <- function(results) {
  
  results %>%
    mutate(method2 = c_method_r2(method, thr.r2)) %>%
    myggplot(aes(par.causal, AUC_mean, fill = method2, color = method2)) +
    geom_hline(yintercept = 0.5, linetype = 2) +
    geom_col(position = position_dodge(), alpha = 0.3) +
    geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
                  position = position_dodge(width = 0.9), color = "black", width = 0.2) + 
    scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
                       oob = scales::rescale_none) +
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
         fill = "Method", color = "Method")
}

barplot_causal_all <- function(results) {
  
  results %>%
    barplot_causal() +
    geom_hline(aes(yintercept = auc_max), linetype = 3, color = "blue",
               data = data.frame(par.h2 = c(0.5, 0.8), auc_max = c(0.84, 0.94))) +
    facet_grid(par.dist ~ par.h2)
}

results2.summary %>%
  filter(par.model == "simple") %>%
  barplot_causal_all()
# ggsave("new-results-all.pdf", scale = 1/90, width = 1235, height = 917)


####

results2.summary2 <- results2 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AltSens_mean = mean(AltSens), AltFDP_mean = mean(AltFDP),
            AltSens_boot = boot(AltSens, n = 1e4), AltFDP_boot = boot(AltFDP, n = 1e4))


barplot_causal2 <- function(results, h2 = 0.8) {
  
  results %>%
    filter(par.h2 == h2) %>%
    mutate(method2 = c_method_r2(method, thr.r2)) %>%
    myggplot(aes(par.causal, mean, fill = method2, color = method2)) +
    geom_col(position = position_dodge(), alpha = 0.3) +
    geom_errorbar(aes(ymin = mean - 2 * boot, ymax = mean + 2 * boot),
                  position = position_dodge(width = 0.9), color = "black", width = 0.2) +
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 Sensitivity",
         fill = "Method", color = "Method")
}

results2.summary2 %>% 
  filter(par.model == "simple") %>%
  gather("metric", "value", c("AltSens_mean", "AltFDP_mean")) %>% 
  barplot_causal2() +
  facet_grid(par.dist ~ metric) +
  scale_y_continuous(limits = c(0, 1))
# ggsave("new-results-all2.pdf", scale = 1/90, width = 1235, height = 917)


results2.summary2 %>% 
  filter(par.model == "simple") %>%
  mutate(id = row_number()) %>%
  gather(key = "metric", value = "value", AltSens_mean:AltFDP_boot) %>%
  separate(metric, c("metric", "stat"), sep = "_") %>%
  spread(stat, value) %>% 
  barplot_causal2() +
  facet_grid(par.dist ~ metric) +
  scale_y_continuous(limits = c(0, 1))


####

results2 %>%
  filter(method == "logit-simple") %>%
  group_by_at(c(vars(starts_with("par")), "alpha")) %>%
  summarise(n = n()) %>%
  print(n = Inf)
