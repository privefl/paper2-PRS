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

# Ttrees
write(c(rep("", 5), pheno), ncolumns = 1, tmp <- tempfile())
writeLines(
  system(sprintf("awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' %s %s",
                 tmp, file.jdb), intern = TRUE),
  tmpfile <- tempfile(fileext = ".jdb")
)
cat(ind.train - 1, file = "ttrees_learn.txt", sep = "\t")
cat(ind.test - 1,  file = "ttrees_val.txt",   sep = "\t")
TTree <- "../../TTree-source/TTree"
system(glue::glue(
  "{TTree}",
  " -j {tmpfile}",
  " -m 3",
  " -b ttrees.bloc",
  " -l ttrees_learn.txt",
  " -v ttrees_val.txt", 
  " -t 1000 -k 1000 -c 5 -n 2000",
  " -x -s"
))
