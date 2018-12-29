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

library(lassosum)
unlink("backingfiles/celiac_test*")
celiac$fam$affection <- y
snp_writeBed(celiac, "backingfiles/celiac_test1.bed",
             ind.row = ind.test.split[[1]])
# celiac$fam$affection[] <- 9999      ## just to be sure that this is not used
# snp_writeBed(celiac, "backingfiles/celiac_test2.bed",
#              ind.row = ind.test.split[[2]])
cor <- p2cor(p = pval, n = length(ind.train), sign = gwas$estim)
cor <- ifelse(is.na(cor), max(abs(cor), na.rm = TRUE) * sign(gwas$estim), cor)
library(doParallel)
registerDoParallel(cl <- makeCluster(nb_cores()))
system.time(
  out <- lassosum.pipeline(cor = cor, chr = CHR, pos = POS, A1 = A2, A2 = A1,
                           test.bfile = "backingfiles/celiac_test1",
                           LDblocks = "EUR.hg19",
                           cluster = cl,
                           exclude.ambiguous = FALSE)
) # 64 sec
stopCluster(cl)

v <- validate(out, covar = U[ind.test.split[[1]], ])
length(ind <- which(v$best.beta != 0))
pred <- big_prodVec(G, -v$best.beta[ind], 
                    ind.row = ind.test.split[[2]],  
                    ind.col = ind)

AUC(pred, y[ind.test.split[[2]]])  ## 83.2 [81.1-85.3]



