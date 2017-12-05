results4 <- list.files("results4", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  format_results()

results4 %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC), N = n()) %>%
  print(n = Inf)

results4 <- list.files("results4", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  format_results() %>%
  cbind(simu = "figures/supp-AUC-chr6") %>%
  bind_rows() %>%
  as_tibble()

results34 <- rbind(
  cbind(simu = "all", results3),
  cbind(simu = "chr6", results4)
) %>%
  as_tibble()

results34 %>%
  myggplot() +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(aes(yintercept = auc_max), linetype = 3, color = "blue",
             data = data.frame(par.h2 = c(0.5, 0.8), auc_max = c(0.84, 0.94))) +
  geom_boxplot(aes(method, AUC, fill = simu, color = simu), alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(par.model ~ par.causal) +
  facet_grid(par.h2 + par.dist ~ par.causal) +
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2) + 
  labs(x = "Method")

ggsave("figures/supp-AUC-chr6.pdf", scale = 1/90, width = 1153, height = 908)

# results4.summary %>%
#   filter(par.h2 == 0.8, method %in% c("logit-simple", "PRS-max")) %>%
#   myggplot(aes(par.causal, AUC_mean, fill = method, color = method)) +
#   geom_hline(yintercept = 0.5, linetype = 2) +
#   geom_hline(yintercept = 0.94, linetype = 3, color = "blue") +
#   geom_bar(stat = "identity", alpha = 0.3, position=position_dodge()) +
#   geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
#                 position=position_dodge(width=0.9), color = "black", width = 0.2) +
#   scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
#                      oob = scales::rescale_none) +
#   labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
#        fill = "Method", color = "Method") +
#   facet_grid(par.dist ~ par.model + par.h2)
