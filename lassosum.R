# devtools::install_github("tshmak/lassosum")
library(bigsnpr)
celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
A1 <- celiac$map$allele1
A2 <- celiac$map$allele2
NCORES <- nb_cores()
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

library(lassosum)
snp_writeBed(celiac, "backingfiles/celiac_test1.bed",
             ind.row = ind.test.split[[1]])
snp_writeBed(celiac, "backingfiles/celiac_test2.bed",
             ind.row = ind.test.split[[2]])
cor <- p2cor(p = pval, n = 12e3, sign = gwas$estim)
cor <- ifelse(is.na(cor), max(abs(cor), na.rm = TRUE) * sign(gwas$estim), cor)
library(doParallel)
registerDoParallel(cl <- makeCluster(nb_cores()))
system.time(
  out <- lassosum.pipeline(cor = cor, chr = CHR, pos = POS, A1 = A1, A2 = A2,
                           test.bfile = "backingfiles/celiac_test1",
                           LDblocks = "EUR.hg19",
                           cluster = cl)
) # 64 sec
stopCluster(cl)

v <- validate(out, covar = U[ind.test.split[[1]], ])
out2 <- subset(out, s = v$best.s, lambda = v$best.lambda)
v2 <- validate(out2, covar = U[ind.test.split[[2]], ],
               test.bfile = "backingfiles/celiac_test2")

AUC(v2$pgs[[1]], v2$pheno - 1)  ## 0.204
AUC(v2$best.pgs, v2$pheno - 1)  ## 0.205



