library(bigsnpr)
celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
NCORES <- nb_cores()
y <- celiac$fam$affection - 1

dim(G)
set.seed(1)
ind.train <- sort(sample(nrow(G), size = 12e3))
ind.test <- setdiff(rows_along(G), ind.train)

snp_writeBed(celiac, "backingfiles/celiac_train.bed", ind.row = ind.train)
snp_writeBed(celiac, "backingfiles/celiac_test.bed",  ind.row = ind.test)

system("BOLT-LMM_v2.3.2/bolt -h")
system(paste(
  "BOLT-LMM_v2.3.2/bolt", 
  "--bfile backingfiles/celiac_train", 
  "--phenoUseFam",
  "--lmm",
  "--LDscoresFile=BOLT-LMM_v2.3.2/tables/LDSCORE.1000G_EUR.tab.gz",
  "--statsFile res.txt",
  "--numThreads=6"
)) ## 55 min

gwas <- bigreadr::fread2("res.txt")
all(gwas$BETA != 0)  # TRUE
hist(pval_lmm <- as.numeric(gwas$P_BOLT_LMM))

preds <- big_prodVec(G, -gwas$BETA, ind.test)
AUC(preds, y[ind.test])   # 0.813
