library(tidyverse)
library(bigsnpr)
library(pkg.paper.PRS)

corr <- readRDS("backingfiles/corr2.rds")
system.time(
  results4 <- list.files("results4", full.names = TRUE) %>%
    read_format_results(corr)
)  # 78 sec for 932 of ~150KB

print(results4, width = Inf)

####

results4 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC, n = 1e4), N = n()) %>%
  print(n = Inf) -> results4.summary

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

barplot_causal_all(results4.summary)
# ggsave("new-results.pdf", scale = 1/90, width = 1235, height = 917)


####

results4.summary2 <- results4 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars(starts_with("par")), "method", "thr.r2")) %>%
  summarise(AltSens_mean = mean(AltSens), AltFDP_mean = mean(AltFDP))


barplot_causal2 <- function(results, h2 = 0.8) {
  
  results %>%
    filter(par.h2 == h2) %>%
    mutate(method2 = c_method_r2(method, thr.r2)) %>%
    myggplot(aes(par.causal, value, fill = method2, color = method2)) +
    geom_col(position = position_dodge(), alpha = 0.3) +
    labs(x = "Causal SNPs (number and location)", y = "Mean of 100 Sensitivity",
         fill = "Method", color = "Method")
}

results4.summary2 %>% 
  gather("metric", "value", c("AltSens_mean", "AltFDP_mean")) %>% 
  barplot_causal2() +
  facet_grid(par.dist ~ metric) +
  scale_y_continuous(limits = c(0, 1))


####

results4 %>%
  filter(method == "logit-simple") %>%
  group_by_at(c(vars(starts_with("par")), "alpha")) %>%
  summarise(n = n()) %>%
  print(n = Inf)
