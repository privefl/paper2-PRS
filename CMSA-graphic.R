library(tidyverse)

col_inner <- scales::hue_pal()(2)
col_outer <- c("#d3d3d3", scales::hue_pal()(3)[2])


plot_remains <- function(K = 5, which_remains = K, size = 20) {
  
  ind <- rep(2, K); ind[which_remains] <- 1
  
  data.frame(
    nest = rep(c("outer", "inner"), c(K + 1, 2)),
    Set = c("test", paste0("train-", 1:K), "test", "train"),
    value = c(20, rep(80 / K, K), 20, 80)
  ) %>%
    ggplot(aes(nest, value, fill = Set, 
               label = paste0(value, "%"))) +
    geom_col(width = 1, color = "black") +
    coord_polar("y", start = 0, direction = -1) +
    scale_fill_manual(values = c(col_outer, col_inner[ind])) +
    theme_void(base_size = size)
}

plot_remains(5)
plot_remains(5, size = 30)


p1 <- plot_remains(5, size = 25) +
  geom_text(position = position_stack(vjust = 0.5), size = rel(4)) +
  theme(legend.key.height = unit(1.5, "line"), 
        legend.key.width  = unit(1.5, "line"),
        legend.position = "left")
print(p1)
ggsave("../thesis-docs/figures/CMSA-explained1.svg", p1,
       scale = 1/90, width = 740, height = 580)

source("CMSA-inside.R")
ggsave("../thesis-docs/figures/CMSA-explained2.svg", p2,
       scale = 1/100, width = 900, height = 650)

cowplot::plot_grid(p1, p2, nrow = 1,
                   labels = LETTERS[1:2], 
                   label_size = 25,
                   scale = 0.95, rel_widths = c(4, 7))
ggsave("figures/CMSA-explained.pdf", scale = 1/100, width = 1300, height = 600)
ggsave("../thesis-docs/figures/CMSA-explained.svg", scale = 1/100, width = 1300, height = 600)
