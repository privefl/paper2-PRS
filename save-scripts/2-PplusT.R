library(bigsnpr)
library(ggplot2)

NCORES <- nb_cores()

# data
celiacUK <- snp_attach("backingfiles/celiacQC_sub4.rds")
G <- celiacUK$genotypes 
n <- nrow(G)
m <- ncol(G)
CHR <- celiacUK$map$chromosome
POS <- celiacUK$map$physical.pos

# parameters
h2 <- 0.8 # heritability
M <- 50
K <- 0.3
ind.train <- sort(sample(n, size = 6000))
ind.test <- setdiff(1:n, ind.train)
model <- "gaussian"
inchr <- c(2, 6)
ind.possible <- snp_indLRLDR(CHR, POS, subset(LD.wiki34, ID == "hild12"))

# simulation 
set <- sample(ind.possible, size = M)
effects <- `if`(model == "gaussian", 
                rnorm(M, sd = sqrt(h2 / M)),
                rmutil::rlaplace(M, s = sqrt(h2 / (2*M))))
y.simu <- scale(G[, set]) %*% effects
y.simu <- y.simu / sd(y.simu) * sqrt(h2)
print(var(y.simu))
y.simu <- y.simu + rnorm(n, sd = sqrt(1 - h2))
pheno <- as.numeric(y.simu > qnorm(1 - K))

# gwas
print(system.time(
  gwas.train <- big_univLogReg(G, pheno[ind.train], 
                               ind.train = ind.train, 
                               covar.train = obj.svd2$u[ind.train, ], 
                               ncores = NCORES)
))
# snp_manhattan(gwas.train, infos.chr = CHR, infos.pos = POS,
#               ind.highlight = set, npoints = 20e3)

# snp_qq(gwas.train) + xlim(1, NA) + 
#   aes(color = as.factor(1:m %in% set)) + 
#   labs(color = "Causal?") +
#   scale_color_manual(values = c("black", "red"))

# clumping
ind.keep <- snp_clumping(G, infos.chr = CHR,
                         ind.row = ind.train,
                         S = abs(gwas.train$score),
                         ncores = NCORES)


# PRS
thrs <- c(0, -log10(5e-8), exp(seq(log(0.1), log(100), length.out = 100)))
lpS <- -predict(gwas.train)
nb.pred <- sapply(thrs, function(thr) sum(lpS[ind.keep2] > thr))
print(system.time(
  prs <- snp_PRS(G, betas.keep = gwas.train$estim[ind.keep],
                 ind.test = ind.test,
                 ind.keep = ind.keep,
                 lpS.keep = lpS[ind.keep], 
                 thr.list = thrs)
))
dim(prs)

# Voir Hmisc::::rcorr.cens
# print(aucs <- apply(prs, 2, AUCBoot, target = pheno[ind.test]))
print(aucs2 <- apply(prs, 2, AUC, target = pheno[ind.test]))
aucs2[c(1:2, which.max(aucs2))]

bigstatsr:::MY_THEME(qplot(nb.pred, aucs["Mean", ], main = "PRS", 
                                 xlab = "# of predictors", 
                                 ylab = "AUC (SD)")) + 
        geom_errorbar(aes(ymin = aucs["Mean", ] - aucs["Sd", ], 
                          ymax = aucs["Mean", ] + aucs["Sd", ])) + 
  scale_x_log10()
