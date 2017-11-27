results_ttrees <- results1 %>%
  filter(method %in% c("logit-simple", "T-Trees")) %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(
    AUC_mean = mean(AUC), 
    AUC_boot = boot(AUC, 1e4, mean), 
    N = n(),
    nb_preds_mean = mean(nb.preds),
    timing_mean = mean(timing)
  ) %>%
  print(n = Inf)


results_ttrees 
rle(results_ttrees$N)

myggplot(results_ttrees, aes(par.causal, AUC_mean, fill = method)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_bar(stat = "identity", alpha = 0.5, position = "dodge") +
  geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
                position=position_dodge(width = 0.9), color = "black", width = 0.2) +
  facet_grid(par.model + par.h2 ~ par.dist) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2))) +
  scale_y_continuous(limits = c(0.5, 1), minor_breaks = 0:20 / 20,
                     # https://stackoverflow.com/a/10365218/6103040
                     oob = scales::rescale_none) + 
  scale_fill_brewer(type = "qual", palette = 2) +
  labs(x = "Causal SNPs (number and location)", y = "Mean of 5 AUCs", fill = "Method")

ggsave("figures/ttrees.pdf", scale = 1/90, width = 560, height = 480)