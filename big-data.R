library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes

# G2 <- big_copy(G, ind.row = rep(rows_along(G), 10), 
#                backingfile = "backingfiles/celiacQC_rep10", save = TRUE)
G2 <- big_attach("backingfiles/celiacQC_rep10.rds")

n <- nrow(G)
ind_train <- sort(sample(n, 12e3))
ind_test <- setdiff(ind_unique, ind_train)

ind_train10 <- as.vector(outer(ind_train, 0:9 * n, '+'))
ind_test10  <- as.vector(outer(ind_test,  0:9 * n, '+'))
y_10 <- rep(celiac$fam$affection - 1, 10)

system.time(
  cmsa.logit <- big_spLogReg(X = G2, y01.train = y_10[ind_train10], 
                             ind.train = ind_train10, 
                             ncores = 6, alphas = 0.5)
)
