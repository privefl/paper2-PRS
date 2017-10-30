
library(tidyverse)

myggplot <- function(..., coeff = 1) {
  bigstatsr:::MY_THEME(ggplot(...), coeff = coeff)
} 

top10 <- 1:110
top20 <- 1:220

results <- list.files("../thesis/paper2-PRS/results1/", full.names = TRUE) %>%
  map_dfr(~readRDS(.x)) %>%
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
  filter(method == "T-Trees") %>%
  pull(timing)

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
  filter(!grepl("PRS", method)) %>%
  plot_results("timing", "Timing") +
  scale_y_continuous(breaks = 0:10 * 2000, minor_breaks = NULL) 
ggsave("timing.pdf", scale = 1/90, width = 1200, height = 900)

plot_results(results, "AUC") +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)
ggsave("AUC.pdf", scale = 1/90, width = 1200, height = 900)

ggsave("perc-cases.pdf", scale = 1/90, width = 1200, height = 900)

plot_results(results, "percCases20", "Percentage of cases in top 20%") +
  scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)

plot_results(results, "nb.preds", "Number of predictors") +
  scale_y_log10(breaks = c(10^(0:7), 3 * 10^(0:7)), minor_breaks = NULL,
                labels = scales::comma_format())
ggsave("nb-preds.pdf", scale = 1/90, width = 1200, height = 900)

ttrees_vs_logit <- filter(results, method %in% c("T-Trees", "logit-simple"))

p_list <- list(
  plot_results(ttrees_vs_logit, "timing", "Timing") +
    scale_y_continuous(breaks = 0:10 * 2000, minor_breaks = NULL),
  plot_results(ttrees_vs_logit, "nb.preds", "Number of predictors") +
    scale_y_log10(breaks = c(10^(0:7), 3 * 10^(0:7)), minor_breaks = NULL,
                  labels = scales::comma_format()),
  plot_results(ttrees_vs_logit, "AUC") +
    scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10),
  plot_results(ttrees_vs_logit, "percCases10", "Percentage of cases in top 10%") +
    scale_y_continuous(breaks = 0:10 / 10, minor_breaks = c(0:9 + 0.5) / 10)
)
lapply(p_list, function(p) p + theme(legend.position = "none")) %>%
  cowplot::plot_grid(plotlist = ., ncol = 2, align = "hv", scale = 0.9,
                     labels = LETTERS[1:4], label_size = 25) %>%
  cowplot::plot_grid(cowplot::get_legend(p_list[[1]]),
                     rel_widths = c(1, 0.15))

