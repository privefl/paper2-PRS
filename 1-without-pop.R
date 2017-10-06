library(bigsnpr)
library(ggplot2)

celiac <- snp_attach("backingfiles/celiacQC.rds")

# Compute PCA
obj.svd <- snp_autoSVD(celiac$genotypes, 
                       celiac$map$chromosome,
                       celiac$map$physical.pos, 
                       ncores = nb_cores())

plot(obj.svd, type = "scores", scores = 1:2) + 
  aes(color = celiac$fam$family.ID) + 
  labs(color = "Population")

# Which are UK controls?
ind.UKcontrols <- which(startsWith(celiac$fam$family.ID, "UK") & 
                          celiac$fam$affection == 1)
obj.svd$u[-ind.UKcontrols, ] <- NA

# Compute Mahalanobis distance on PCs
distM <- rep(NA, nrow(celiac$genotypes))
distM[ind.UKcontrols] <- robust::covRob(obj.svd$u[ind.UKcontrols, ])$dist
THR <- 30
# Plot UK controls, coloring with Mahalanobis distance
plotlist <- lapply(list(1:2, 3:4, 5:6, 7:8), function(scores) {
  plot(obj.svd, type = "scores", scores = scores, coeff = 0.7) + 
    aes(color = (distM > THR)) + 
    labs(color = "Outlier?")
}) 
cowplot::plot_grid(plotlist = plotlist, ncol = 2, align = "hv")

# New dataset with only non-outliers UK controls
length(ind.keep <- which(distM < THR))
celiac2 <- snp_attach(subset(celiac, ind.row = ind.keep))
obj.svd2 <- snp_autoSVD(celiac2$genotypes, 
                        celiac2$map$chromosome,
                        celiac2$map$physical.pos,
                        ncores = nb_cores())
plot(obj.svd2)
plot(obj.svd2, type = "scores", scores = 1:2) 
plot(obj.svd2, type = "scores", scores = 3:4)
