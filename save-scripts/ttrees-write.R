library(bigsnpr)
library(ggplot2)

NCORES <- nb_cores()

# data
celiac2 <- snp_attach("backingfiles/celiacQC_sub4.rds")
G <- celiac2$genotypes 
n <- nrow(G)
m <- ncol(G)
CHR <- celiac2$map$chromosome
POS <- celiac2$map$physical.pos

# jdb
file.jdb <- "ttrees.jdb"
cat(rep(";", 5), sep = "\n", file = file.jdb)
big_apply(G, a.FUN = function(X, ind) {
  mat <- cbind(sample(0:1, length(ind), replace = TRUE), X[ind, ])
  write.table(mat, file = file.jdb, append = TRUE, quote = FALSE,
              row.names = paste0("ind_", ind), col.names = FALSE)
  NULL
}, a.combine = "c", ind = rows_along(G), block.size = 1e3)
# bloc
bigsnpr:::write.table2(bigstatsr:::CutBySize(ncol(G), block.size = 10)[, 1:2] - 1,
                       file = "ttrees.bloc")

