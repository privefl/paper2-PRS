library(bigsnpr)
obj.svd <- readRDS("backingfiles/PCA.rds")
plot(obj.svd, type = "scores")
# Get population
celiac <- snp_attach("backingfiles/celiacQC.rds")

library(ggplot2)
p1 <- plot(obj.svd, type = "scores", scores = 1:2, coeff = 1.2) + 
  aes(color = celiac$fam$family.ID) +
  labs(color = "Population") +
  guides(colour = guide_legend(override.aes = list(size = 4)))
p2 <- plot(obj.svd, type = "scores", scores = 3:4, coeff = 1.2) + 
  aes(color = celiac$fam$family.ID) +
  labs(color = "Population")

library(magrittr)
library(cowplot)
plot_grid(
  plotlist = lapply(list(p1, p2), function(p) {
    p + theme(legend.position = "none") + ggtitle(NULL)
  })
) %>%
  plot_grid(get_legend(p1), rel_widths = c(1, 0.15), scale = 0.95)
  
ggsave("../thesis-docs/figures/PC1-4.png", scale = 1/90, width = 1270, height = 610)
