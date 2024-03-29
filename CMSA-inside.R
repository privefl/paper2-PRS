# system.time(
#   logit <- big_spLogReg(X = G, y01.train = y[ind.train],
#                         ind.train = ind.train,
#                         covar.train = obj.svd$u[ind.train, ],
#                         ncores = NCORES,
#                         nlam.min = 150,
#                         return.all = TRUE)
# )

myggplot <- function(..., coeff = 1) {
  bigstatsr:::MY_THEME(ggplot(...), coeff = coeff) +
    theme(strip.text.x = element_text(size = rel(2)),
          strip.text.y = element_text(size = rel(2)))
} 

# tmp <- logit[[2]]
# saveRDS(tmp, "one-logit-CMSA.rds")
tmp <- readRDS("one-logit-CMSA.rds")
ind.min <- which.min(tmp$loss.val)
ind <- seq(2, length(tmp$lambda), by = 2)

library(tidyverse)

p2 <- bind_rows(
  data_frame(lambda = tmp$lambda[ind], loss = tmp$loss.val[ind], 
             Dataset = "Remaining fold"),
  data_frame(lambda = tmp$lambda[ind], loss = tmp$loss[ind],
             Dataset = "Training folds")
) %>%
  myggplot(coeff = 1.2) +
  geom_vline(xintercept = c(tmp$lambda[ind.min + c(0, 10)]), 
             linetype = c(2, 3), size = 1) +
  geom_point(aes(lambda, loss, color = Dataset)) +
  scale_x_log10(limits = c(0.006, NA)) +
  labs(x = "Lambda\n(path: right to left)", y = "Loss") +
  theme(legend.position = c(0.75, 0.25)) + 
  guides(colour = guide_legend(override.aes = list(size = 3)))

print(p2)

