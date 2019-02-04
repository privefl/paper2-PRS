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
pval[pval == 0] <- .Machine$double.xmin

library(dplyr)
celiac$map %>%
  select(hg19chrc = chromosome, snpid = marker.ID,
         a1 = allele1, a2 = allele2, bp = physical.pos) %>%
  mutate(hg19chrc = paste0("chr", hg19chrc), or = exp(gwas$estim), p = pval) %>%
  # filter(hg19chrc %in% c("chr1", "chr2")) %>%
  bigreadr::fwrite2("SUM_STATS_FILE.txt", sep = "\t")

readLines("SUM_STATS_FILE.txt", n = 5)
# [1] "hg19chrc\tsnpid\ta1\ta2\tbp\tor\tp"                                 
# [2] "chr1\trs3934834\tA\tG\t995669\t1.03012433299041\t0.449674955932202" 
# [3] "chr1\trs3737728\tA\tG\t1011278\t1.02982046003153\t0.363400424211927"
# [4] "chr1\trs6687776\tA\tG\t1020428\t1.08042108970176\t0.048125118630458"
# [5] "chr1\trs9651273\tA\tG\t1021403\t1.0347803708118\t0.291137612222826"

unlink(paste0("backingfiles/celiac_test1", c(".bed", ".bim", ".fam")))
snp_writeBed(celiac, "backingfiles/celiac_test1.bed", ind.row = ind.test.split[[1]])

reticulate::use_python("/home/privef/anaconda3/bin/python3")
reticulate::py_config()
stopifnot(system("python3 --version", intern = TRUE) == "Python 3.7.0")
ldpred <- "../ldpred/LDpred.py"
unlink("OUT_COORD_FILE.hdf5")
system(glue::glue(
  "python3 {ldpred} coord", 
  " --gf backingfiles/celiac_test1",
  " --ssf SUM_STATS_FILE.txt --ssf-format BASIC",
  " --N {length(ind.train)}",
  " --out OUT_COORD_FILE.hdf5"
))

system(glue::glue(
  "ldpred", 
  " --coord OUT_COORD_FILE.hdf5",
  " --ld_radius {round(ncol(G) / 3000)}",
  " --PS 0.01,0.001,0.0001",
  " --N {length(ind.train)}",
  " --out OUT_COORD_FILE"
))

res <- bigreadr::fread2("OUT_COORD_FILE.txt_LDpred_p1.0000e-04.txt")
plot(ldpred_beta ~ raw_beta, data = res, pch = 20)
