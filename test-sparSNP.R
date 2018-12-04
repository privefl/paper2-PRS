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
snp_writeBed(celiac, "backingfiles/celiac_test.bed", ind.row = ind.test)

Sys.setenv(PATH = paste0("~/Bureau/SparSNP:", Sys.getenv("PATH")))
path <- setwd("~/Bureau/SparSNP")
unlink("discovery", recursive = TRUE)
system.time(
  system(sprintf("./crossval.sh %s sqrhinge",
                 file.path(path, "backingfiles/celiac_train")))
)


       