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

library(dplyr)
celiac$map %>%
  select(hg19chrc = chromosome, snpid = marker.ID,
         a1 = allele1, a2 = allele2, bp = physical.pos) %>%
  mutate(or = exp(gwas$estim), p = pval) %>%
  bigreadr::fwrite2("SUM_STATS_FILE.txt", sep = "\t")
  head()
# hg19chrc    snpid    a1    a2    bp    or    p
# chr1    rs4951859    C    G    729679    0.97853    0.2083
# chr1    rs142557973    T    C    731718    1.01949    0.3298


system(glue::glue(
  "coord", 
  " --gf=backingfiles/celiac_test1",
  " --ssf=SUM_STATS_FILE",
  " --SSF_FORMAT=BASIC",
  " --N={length(ind.train)}",
  " --out=OUT_COORD_FILE"
))
