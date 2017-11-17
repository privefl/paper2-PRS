# Put all results in a single tibble
results3 <- list.files("results3", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  as_tibble() %>%
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2]))
  )

results3 %>%
  filter(par.h2 == 0.8) %>%
  plot_results("AUC") +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)
