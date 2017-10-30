
library(tidyverse)

myggplot <- function(..., coeff = 1) {
  bigstatsr:::MY_THEME(ggplot(...), coeff = coeff)
} 

top10 <- 1:110
top20 <- 1:220

results <- list.files("../thesis/paper2-PRS/results1/", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
  filter(method != "T-Trees") %>%
  as_tibble() %>%
  mutate(
    par.causal = factor(map_chr(par.causal, ~paste(.x[1], .x[2], sep = " in ")),
                        levels = c("30 in HLA", paste(3 * 10^(1:3), "in all"))),
    AUC = map_dbl(eval, ~bigstatsr::AUC(.x[, 1], .x[, 2])),
    percCases10 = map_dbl(eval, ~mean(.x[order(.x[, 1], decreasing = TRUE)[top10], 2])),
    percCases20 = map_dbl(eval, ~mean(.x[order(.x[, 1], decreasing = TRUE)[top20], 2]))
  )

print(results, n = Inf)


results %>%
  group_by_at(c(vars(starts_with("par")), "method")) %>%
  summarise_at(c("timing", "nb.preds", "AUC", "percCases10"), mean)

plot_results <- function(results, y, ylab = y) {
  
  dist <- "Distribution\nof effects"
  
  myggplot(results) +
    geom_boxplot(aes_string("method", y, color = "par.dist", 
                            fill = "par.dist"), alpha = 0.3) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    facet_grid(par.model ~ par.causal) +
    theme(strip.text.x = element_text(size = rel(2)),
          strip.text.y = element_text(size = rel(2))) +
    labs(x = "Method", y = ylab, fill = dist, color = dist)
}


results %>%
  filter(par.h2 == 0.8) %>%
  plot_results("AUC") +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)
ggsave("../thesis/paper2-PRS/figures/AUC.pdf", 
       scale = 1/90, width = 1200, height = 900)

ggsave("perc-cases.pdf", scale = 1/90, width = 1200, height = 900)

results %>%
  filter(par.h2 == 0.8) %>%
  plot_results("percCases20", "Percentage of cases in top 20%") +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)

results %>%
  filter(par.h2 == 0.8) %>%
  plot_results("nb.preds", "Number of predictors") +
  scale_y_log10(breaks = c(10^(0:7), 3 * 10^(0:7)), minor_breaks = NULL,
                labels = scales::comma_format())
ggsave("nb-preds.pdf", scale = 1/90, width = 1200, height = 900)

results %>%
  filter(grepl("logit", method)) %>%
  group_by(method, par.causal) %>%
  summarise_at("nb.preds", list("Min" = min, "Max" = max))


#### Correlation between 


ggsave("../thesis/paper2-PRS/figures/AUC-corr.pdf", 
       scale = 1/90, width = 1300, height = 950)

summary(lm(AUC ~ percCases10, data = results))
summary(lm(AUC ~ percCases20, data = results))

results %>%
  group_by_at(c(vars(starts_with("par")), "num.simu")) %>%
  mutate(AUC_rel = AUC / AUC[method == "logit-simple"]) %>%
  # select(-eval, -nb.preds, -percCases10, -percCases20) %>%
  plot_results(y = "AUC_rel", ylab = "Relative AUC / 'logit-simple'")

