# devtools::install_github("amkusmec/FarmCPUpp")
library(FarmCPUpp)
library(bigmemory)
library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
# ind.chr6 <- which(celiac$map$chromosome == 6)
# 
# G2 <- filebacked.big.matrix(nrow(G), ncol(G),
#                             backingfile = "test_farmcpu")
# 
# big_apply(G, a.FUN = function(X, ind) {
#   G2[, ind] <- X[, ind]
#   NULL
# }, a.combine = 'c')
# G2[, ncol(G2)]

G2 <- attach.big.matrix("test_farmcpu.desc")

system.time(
  myResults2 <- farmcpu(
    Y = celiac$fam[c(1, 6)], 
    GD = G2, 
    GM = celiac$map[c(2, 1, 4)],
    CV = readRDS("backingfiles/PCA.rds")$u,
    ncores.glm = nb_cores(),
    ncores.reml = nb_cores()
  )
) 
# 69 / 42 for chr 6
# 1880 / 414 for all chrs
# With PCs: 2255 / 479

plot(-log10(myResults$affection$GWAS$p.value),
     -log10(myResults2$affection$GWAS$p.value),
     pch = 20, cex = 0.5)

str(myResults2)
plot(-log10(myResults2$affection$GWAS$p.value), pch = 20, cex = 0.5, ylim = c(0, 30))
plot(-log10(ppoints(length(myResults2$affection$GWAS$p.value))),
     sort(-log10(myResults2$affection$GWAS$p.value), decreasing = TRUE),
     pch = 20, cex = 0.5, ylim = c(0, 30))
abline(0, 1, col = "red")
