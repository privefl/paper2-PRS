# system.time(
#   logit <- big_spLogReg(X = G, y01.train = y[ind.train], 
#                         ind.train = ind.train, 
#                         covar.train = obj.svd$u[ind.train, ],
#                         ncores = NCORES, 
#                         nlam.min = 150,
#                         return.all = TRUE)
# )


tmp <- logit[[2]]

library(tidyverse)

bind_rows(
  data_frame(lambda = tmp$lambda[-1], loss = tmp$loss.val[-1], 
             Dataset = "Remaining fold"),
  data_frame(lambda = tmp$lambda[-1], loss = tmp$loss[-1],
             Dataset = "Training folds")
) %>%
  myggplot(coeff = 1.2) +
  geom_point(aes(lambda, loss, color = Dataset)) +
  scale_x_log10(limits = c(0.006, NA)) +
  labs(x = "Lambda (right to left)", y = "Loss for the remaining fold") +
  theme(legend.position = c(0.75, 0.25))
