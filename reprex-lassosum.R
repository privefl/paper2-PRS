# devtools::install_github("privefl/bigsnpr")
library(bigsnpr)

set.seed(1)
N <- 1000; M <- 1e4; m <- 20
snp <- snp_fake(N, M)
snp$genotypes[] <- sample(0:2, size = N * M, replace = TRUE)
snp$map$physical.pos <- seq(2, 2e8, length.out = M)
set <- sample(M, m)
beta <- rnorm(m) / sqrt(m)
s <- drop(scale(snp$genotypes[, set]) %*% beta)
y <- ((s + rnorm(m)) > quantile(s, 0.7)) + 0
AUC(s, y)
snp$fam$affection <- y + 1

bed <- snp_writeBed(snp, tempfile(fileext = ".bed"))
rds <- snp_readBed(bed, sub("\\.bed$", "", bed))

# devtools::install_github("tshmak/lassosum")
snp <- snp_attach(rds)
G <- snp$genotypes
CHR <- snp$map$chromosome
POS <- snp$map$physical.pos
A1 <- snp$map$allele1
A2 <- snp$map$allele2

set.seed(1)
ind.train <- sort(sample(nrow(G), size = N * 0.6))
ind.test <- setdiff(rows_along(G), ind.train)
ind.test.split <- split(ind.test, sample(1:2, length(ind.test), TRUE))

gwas <- big_univLogReg(G, y[ind.train], ind.train)
plot(gwas, type = "Manhattan")

library(lassosum)
bed2 <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                     ind.row = ind.test.split[[1]])
bed3 <- snp_writeBed(snp, tempfile(fileext = ".bed"),
                     ind.row = ind.test.split[[2]])
cor <- p2cor(p = predict(gwas, log10 = FALSE), n = length(ind.train), 
             sign = gwas$estim)
out <- lassosum.pipeline(cor = cor, chr = CHR, pos = POS, A1 = A2, A2 = A1,
                         test.bfile = sub("\\.bed$", "", bed2),
                         LDblocks = "EUR.hg19")

v <- validate(out)
out2 <- subset(out, s = v$best.s, lambda = v$best.lambda)
v2 <- validate(out2, test.bfile = sub("\\.bed$", "", bed3))

v2$best.validation.result        ## -0.07595338
AUC(v2$best.pgs, v2$pheno - 1)   ## 0.4441773


