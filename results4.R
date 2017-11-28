results4 <- list.files("results4", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2]))
  ) %>%
  cbind(simu = "chr6") %>%
  bind_rows(cbind(select(results3, -percCases10, -percCases20), simu = "all")) %>%
  as_tibble()

results4 %>%
  filter(par.h2 == 0.8) %>%
  myggplot() +
  geom_boxplot(aes(method, AUC, fill = simu, color = simu), alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(par.model ~ par.causal) +
  theme(strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2))) +
  scale_colour_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  geom_hline(yintercept = 0.94, color = "blue", linetype = 2) +
  facet_grid(par.dist ~ par.model + par.h2 + par.causal)


results4 %>%
  filter(par.h2 == 0.8, method %in% c("logit-simple", "PRS-max"), simu == "chr6") %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC, 1e4, mean)) %>%
  myggplot(aes(par.causal, AUC_mean, fill = method, color = method)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_hline(yintercept = 0.94, linetype = 3, color = "blue") +
  geom_bar(stat = "identity", alpha = 0.3, position=position_dodge()) +
  geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, ymax = AUC_mean + 2 * AUC_boot),
                position=position_dodge(width=0.9), color = "black", width = 0.2) +
  scale_y_continuous(limits = c(0.5, NA), minor_breaks = 0:20 / 20,
                     oob = scales::rescale_none) +
  scale_fill_brewer(type = "qual", palette = 2) +
  scale_color_brewer(type = "qual", palette = 2) +
  labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
       fill = "Method", color = "Method") +
  facet_grid(par.dist ~ par.model + par.h2) + 
  theme(strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2)))
