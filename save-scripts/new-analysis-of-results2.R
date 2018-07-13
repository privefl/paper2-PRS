library(tidyverse)
library(bigsnpr)
library(pkg.paper.PRS)

corr <- readRDS("backingfiles/corr2.rds")
system.time(
  results5 <- list.files("results5", full.names = TRUE) %>%
    read_format_results(corr)
)  # 62 sec for 647 of ~1Mb

results5 %>%
  filter(method %in% c("PRS-max", "logit-simple")) %>%
  group_by_at(c(vars("n.train", starts_with("par")), "method", "thr.r2")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC, n = 1e4), N = n()) %>%
  print(n = Inf) -> results5.summary

results5.summary %>%
  mutate(method2 = c_method_r2(method, thr.r2)) %>%
  ggplot(aes(n.train, AUC_mean, color = method2, linetype = method2)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(aes(yintercept = auc_max), linetype = 3, color = "blue",
             data = data.frame(par.h2 = c(0.5, 0.8), auc_max = c(0.84, 0.94))) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, 
                    ymax = AUC_mean + 2 * AUC_boot), 
                size = 1.5, width = 0) +
  facet_grid(par.h2 ~ par.dist) +
  scale_x_continuous(breaks = 1:6 * 1000, minor_breaks = NULL) + 
  scale_y_continuous(minor_breaks = 0:20 / 20) +
  # scale_color_manual(values = methods.color) +
  labs(x = "Size of the training set", y = "Mean of 100 AUCs", 
       color = "Method", linetype = "Method") +
  theme(legend.key.width = unit(3, "line")) +
  theme_bigstatsr()


