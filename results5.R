results5 <- list.files("results5", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  as_tibble() %>%
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2]))
  )

pryr::object_size(results5)

results5_3 <- results3 %>%
  select(-eval) %>%
  filter(par.model == "simple", par.causal == "300 in all") %>%
  cbind(n.train = 6000, .) %>%
  bind_rows(select(results5, -eval)) %>%
  as_tibble()

results5_3 %>%
  group_by_at(c(vars(starts_with("par")), "method", "n.train")) %>%
  summarise(AUC_mean = mean(AUC), AUC_boot = boot(AUC, 1e4, mean), N = n()) %>%
  print(n = Inf)

tmp <- .Last.value

myggplot(tmp, aes(n.train, AUC_mean, color = method)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  geom_point(size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(ymin = AUC_mean - 2 * AUC_boot, 
                    ymax = AUC_mean + 2 * AUC_boot), 
                size = 1.5, width = 0) +
  facet_grid(par.model + par.h2 ~ par.causal + par.dist) +
  theme(strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2))) +
  scale_x_continuous(breaks = 1:6 * 1000, minor_breaks = NULL) + 
  scale_y_continuous(minor_breaks = 0:20 / 20) +
  scale_fill_brewer(type = "qual", palette = 2) +
  labs(x = "Causal SNPs (number and location)", y = "Mean of 100 AUCs",
       fill = "Method")
