res_auc <- data.frame(AUC = c(0.587, 0.587), Method = c("C+T-max", "PLR"))
res_cor <- data.frame(cor = c(0.555, 0.655), Method = c("C+T-max", "PLR"))
library(ggplot2)
p1 <- ggplot() +
  geom_col(aes(Method, AUC, fill = Method), data = res_auc, position = "dodge",
           width = 0.9, alpha = 0.5, color = "black", size = 1) +
  geom_hline(yintercept = 0.5, linetype = 2) +
  bigstatsr::theme_bigstatsr() +
  scale_y_continuous(limits = c(0.5, 0.6), breaks = seq(0.5, 0.6, 0.02),
                     minor_breaks = seq(0.5, 0.6, 0.01), oob = scales::rescale_none) +
  ggtitle("PRS for Breast Cancer")


p2 <- ggplot() +
  geom_col(aes(Method, cor, fill = Method), data = res_cor, position = "dodge",
           width = 0.9, alpha = 0.5, color = "black", size = 1) +
  bigstatsr::theme_bigstatsr() +
  scale_y_continuous(breaks = seq(0, 1, 0.05), minor_breaks = seq(0, 1, 0.01)) +
  labs(y = "Correlation (within each sex)") +
  ggtitle("PRS for Height")

library(cowplot)
plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  scale = 0.95
)

