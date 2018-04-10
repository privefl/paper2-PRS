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
    ggplot(aes(nest, value, fill = Set)) +
    geom_col(width = 1, color = "black") +
    coord_polar("y", start = 0, direction = -1) +
    scale_fill_manual(values = c(col_outer, col_inner[ind])) +
    theme_void(base_size = size)
}

plot_remains(5)
plot_remains(5, size = 25)

