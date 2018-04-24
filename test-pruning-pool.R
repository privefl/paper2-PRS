library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes
CHR <- celiac$map$chromosome
POS <- celiac$map$physical.pos
obj.svd <- readRDS("backingfiles/PCA.rds")
# plot(obj.svd, type = "loadings", loadings = 1:10, coeff = 0.5)

files <- list.files("backingfiles/others/", "^cluster",
           full.names = TRUE)
pop <- snp_getSampleInfos(celiac, files)[[1]]
pop[pop == 4] <- 3
library(ggplot2)
plot(obj.svd, type = "scores", scores = 3:4) +
  aes(color = as.factor(pop))

ind.pop <- split(rows_along(G), pop)

afs <- sapply(seq_along(ind.pop), function(i) {
  big_scale()(G, ind.row = ind.pop[[i]])$center / 2
})

# bigsnpr:::write.table2(t(afs), "backingfiles/celiac-pool.txt")

library(pcadapt)
pool_mat <- read.pcadapt(t(afs), type = "pool")
dim(pool_mat)
class(pool_mat)

obj.pool <- pcadapt(pool_mat)
plot(obj.pool)
scree_plot(obj.pool, 4)
plot(obj.pool$loadings[, 1], pch = 20)
plot(obj.pool$loadings[, 2], pch = 20)
plot(obj.pool$loadings[, 3], pch = 20)
plot(obj.pool$loadings[, 4], pch = 20)

plink <- download_plink()
file.base <- "../Dubois2010_data/FinnuncorrNLITUK1UK3hap300_QC_norel"
tmp <- tempfile()
system(glue::glue("{plink} --bfile {file.base} --out {tmp}", 
                  " --indep-pairwise 500 1 0.05"))
to_keep <- read.table(paste0(tmp, ".prune.in"), stringsAsFactors = FALSE)[[1]]
ind.keep <- match(to_keep, celiac$map$marker.ID)

pool_mat.sub <- read.pcadapt(t(afs[ind.keep, ]), type = "pool")
pcadapt.sub <- pcadapt(pool_mat.sub)
plot(pcadapt.sub, chr.info = celiac$map$chromosome[ind.keep])
scree_plot(pcadapt.sub, 5)
plot(pcadapt.sub$loadings[, 1], pch = 20)
plot(pcadapt.sub$loadings[, 2], pch = 20)
plot(pcadapt.sub$loadings[, 3], pch = 20)

obj.pcadapt <- snp_pcadapt(G, obj.svd$u)
snp_manhattan(obj.pcadapt, CHR, POS, npoints = 50e3)

tmat <- scale(t(afs))
svd <- svd(tmat[, ind.keep])
w <- crossprod(tmat, sweep(svd$u, 2, svd$d, '/')[, 1:3])
plot(w[ind.keep, ], svd$v[, 1:3], pch = 20)
abline(0, 1, col = "red")

res <- pcadapt:::get_statistics(w, method = "mahalanobis", 
                                pass = rows_along(w))
plot(-log10(res$pvalues), pch = 20)

obj.pcadapt.sub <- snp_pcadapt(G, obj.svd$u, ind.col = ind.keep)
snp_manhattan(obj.pcadapt.sub, CHR[ind.keep], POS[ind.keep])
