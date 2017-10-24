
library(tidyverse)

myggplot <- function(..., coeff = 1) {
  bigstatsr:::MY_THEME(ggplot(...), coeff = coeff)
} 

results <- list.files("../thesis/paper2-PRS/results1/", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  as_tibble() %>%
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2]))
  )
  
print(results, n = Inf)

results %>%
  filter(method == "T-Trees") %>%
  pull(timing)

results %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise_at(c("timing", "nb.preds", "AUC"), mean)

results %>%
  filter(!grepl("PRS", method)) %>%
  myggplot() +
  geom_boxplot(aes(method, timing, color = par.dist, fill = par.dist), alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(par.model ~ par.causal) +
  scale_y_continuous(breaks = 0:10 * 2000, minor_breaks = NULL) + 
  theme(strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2))) 

results %>%
  myggplot() +
  geom_boxplot(aes(method, AUC, color = par.dist, fill = par.dist), alpha = 0.3) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  facet_grid(par.model ~ par.causal) +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = NULL) + 
  theme(strip.text.x = element_text(size = rel(2)),
        strip.text.y = element_text(size = rel(2))) 
