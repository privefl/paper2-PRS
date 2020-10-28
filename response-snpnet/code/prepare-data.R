##### Prepare bed file (merge + QC) ####

# file.symlink("~/NCRR-PRS/faststorage/UKBB/", ".")

write(sapply(1:22, function(chr) {
  paste0("UKBB/bed/",
         c(paste0("ukb_cal_chr", chr, "_v2.bed"),
           paste0("ukb_snp_chr", chr, "_v2.bim"),
           paste0("ukb58024_cal_chr", chr, "_v2_s488264.fam")))
}), tmp <- tempfile(), ncolumns = 3)

library(bigsnpr)
snp_plinkQC(
  plink.path = download_plink("tmp-data"),
  prefix.in = tmp,
  file.type = "--merge-list",
  prefix.out = "data/ukbb",
  geno = 0.01,
  mind = 0.1,
  maf = 0.01,
  autosome.only = TRUE,
  extra.options = "--memory 100000"
)
# 504,139 variants and 488,371 people pass filters and QC.


#### Prepare subset, phenotypes and covariates ####

library(dplyr)
library(bigreadr)
csv <- "UKBB/ukb41181.csv"
df0 <- fread2(
  csv,
  select = c("eid", "50-0.0", "21001-0.0", "21022-0.0", "22001-0.0", "22006-0.0",
             "22020-0.0", paste0("22009-0.", 1:16)),
  col.names = c("eid", "height", "BMI", "age", "sex", "white_british",
                "used_in_pca", paste0("PC", 1:16))
) %>%
  filter(!is.na(eid), white_british == 1, used_in_pca == 1) %>%
  select(-white_british, -used_in_pca)
# 337,475 individuals remaining

# Sets: #5 will be test (20%), 3 of them will be training (60%)
#       and the last will be validation (20%)
df0$set <- sample(rep_len(1:5, nrow(df0)))

# Read from bed to bk/rds
system.time(
  snp_readBed2("data/ukbb.bed",
               ind.row = match(df0$eid, fread2("data/ukbb.fam")[[2]]),
               ncores = 16)
) # 6 min

# Impute with mean (+ save with new code.256 CODE_DOSAGE)
ukb <- snp_attach("data/ukbb.rds")
system.time(
  G2 <- snp_fastImputeSimple(ukb$genotypes, method = "mean2", ncores = 16)
) # 5 min
ukb$genotypes <- G2
snp_save(ukb)


#### Other phenotypes ####

ind_sub <- match(df0$eid, fread2(csv, select = "eid")$eid)

df_illness <- fread2(csv, select = c(paste0("20002-0.", 0:28),
                                     paste0("20002-1.", 0:28),
                                     paste0("20002-2.", 0:28)))[ind_sub, ]
df_ICD10 <- fread2(csv, select = c(paste0("40001-", 0:2, ".0"),
                                   paste0("40002-0.", 0:13),
                                   paste0("40002-1.", 0:13),
                                   paste0("40002-2.", 0:13),
                                   paste0("41202-0.", 0:379),
                                   paste0("41204-0.", 0:434)))[ind_sub, ]

# Asthma
ind_asthma <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x == 1111)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 3) == "J45"))
))))
ind_respi <- sort(unique(unlist(c(
  lapply(df_illness, function(x) which(x %in% 1111:1125)),
  lapply(df_ICD10,   function(x) which(substr(x, 1, 1) == "J"))
))))
y <- rep(0, nrow(df_illness)); y[ind_respi] <- NA; y[ind_asthma] <- 1
table(y, exclude = NULL)
#        0      1   <NA>
#   256752  45256  35467
df0$asthma <- y

# High cholesterol
y <- (rowSums(df_illness == 1473, na.rm = TRUE) > 0) + 0
table(y, exclude = NULL)
#      0      1
# 292256  45219
df0$high_cholesterol <- y

saveRDS(df0, "data/pheno_covar.rds")
