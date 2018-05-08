N <- 1000

test <- replicate(1000, {
  f1 <- runif(1)
  f2 <- runif(1)
  g1 <- rbinom(N, size = 2, prob = f1)
  g1.1 <- g1 > 0.5
  g1.2 <- g1 > 1.5
  g2 <- rbinom(N, size = 2, prob = f2)
  g2.1 <- g2 > 0.5
  g2.2 <- g2 > 1.5
  
  g12 <- scale(cbind(g1, g1.1, g1.2, g2, g2.1, g2.2))
  
  g <- scale(g1) * scale(g2)
  
  g12.sub <- g12[, attr(g12, "scaled:scale") > 1e-4, drop = FALSE]
  mod <- glmnet::glmnet(g12.sub, g, alpha = 0.01)
  # coef(mod)
  
  pred <- predict(mod, g12.sub)
  corr <- apply(pred, 2, cor, g)
  # plot(pred[, which.max(corr)], g)
  
  c(f1, f2, max(corr[-1]))
})

plot(test[1, ], test[3, ], pch = 20)
plot(abs(test[1, ] - test[2, ])^2, test[3, ], pch = 20)
rgl::plot3d(test[1, ], test[2, ], test[3, ])

library(ggplot2)
qplot(test[1, ], test[2, ], col = test[3, ]) + 
  viridis::scale_color_viridis() + 
  coord_equal()
