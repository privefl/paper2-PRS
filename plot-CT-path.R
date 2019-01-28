# devtools::install_github("tshmak/lassosum")
library(bigsnpr)
celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
A1 <- celiac$map$allele1
A2 <- celiac$map$allele2
y <- celiac$fam$affection - 1

dim(G)
set.seed(1)
ind.train <- sort(sample(nrow(G), size = 12e3))
ind.test <- setdiff(rows_along(G), ind.train)
ind.test.split <- split(ind.test, sample(1:2, length(ind.test), TRUE))

U <- predict(readRDS("backingfiles/PCA.rds"))

system.time(
  gwas <- big_univLogReg(G, y[ind.train], ind.train,
                         covar.train = U[ind.train, ],
                         ncores = nb_cores())
) # 140 sec
hist(pval <- predict(gwas, log10 = FALSE))
lpval <- -predict(snp_gc(gwas))

ind.keep <- snp_clumping(G, CHR, ind.train, S = lpval, ncores = nb_cores())
thrs <- c(0, -log10(5e-8), exp(seq(log(0.1), log(500), length.out = 100)))
scores <- snp_PRS(G, gwas$estim[ind.keep], ind.test, ind.keep,
                  lpS.keep = lpval[ind.keep], thr.list = thrs)
aucs <- apply(scores, 2, AUC, y[ind.test])

library(ggplot2)
ggplot(data.frame(thr = thrs, auc = aucs), aes(thr, auc)) + 
  geom_point() + 
  geom_line() + 
  labs(x = "-log10(p-value) threshold", y = "AUC of C+T model") + 
  bigstatsr::theme_bigstatsr() + 
  scale_x_log10()

AUCBoot(scores[, which.max(aucs)], y[ind.test])
