library(bigsnpr)

# snp_readBed("../2017_project_5/data/hapmap3_qc.bed")
hapmap3 <- snp_attach("../2017_project_5/data/hapmap3_qc.rds")
G <- hapmap3$genotypes
CHR <- hapmap3$map$chromosome
POS <- hapmap3$map$physical.pos

counts <- big_counts(G)
counts[, 1:4]
dim(G)

snp_fastImpute(G, CHR, ncores = nb_cores())

rds <- "../2017_project_5/data/hapmap3_qc-infos-impute.rds"

tmp <- big_attach(rds)
ind.col <- which(!is.na(tmp[1, ]))

G <- G$copy(code = bigsnpr:::CODE_IMPUTE_PRED)
obj.svd <- snp_autoSVD(G, CHR, POS, ind.col = ind.col, k = 20,
                       ncores = nb_cores())
plot(obj.svd) + 
  ggplot2::scale_y_log10()
plot(obj.svd, type = "loadings", loadings = 1:6, coeff = 0.5)
plot(obj.svd, type = "scores", scores = 9:10)

test <- varimax(obj.svd$v) 
u <- obj.svd$u
u_rot <- u %*% test$rotmat

# obj.pcadapt <- snp_pcadapt(G, u, ind.col = ind.col)
# snp_qq(obj.pcadapt) + ggplot2::xlim(1, NA)
# snp_manhattan(snp_gc(obj.pcadapt), 
#               infos.chr = CHR[ind.col],
#               infos.pos = POS[ind.col],
#               npoints = 50e3) + 
#   ggplot2::geom_hline(yintercept = -log10(5e-8), color = "red")

infos <- data.table::fread("ftp://ftp.ncbi.nlm.nih.gov/hapmap/genotypes/hapmap3/plink_format/draft_1/relationships_w_pops_051208.txt", 
                  data.table = FALSE)

pop <- as.factor(infos$population[match(hapmap3$fam$sample.ID, infos$IID)])

plot(u[, 1:2], col = pop, pch = 20)
plot(u_rot[, 1:2], col = pop, pch = 20)
plot(u_rot[, 3:4], col = pop, pch = 20)
plot(u_rot[, 5:6], col = pop, pch = 20)
plot(u_rot[, 7:8], col = pop, pch = 20)
plot(u_rot[, 9:10], col = pop, pch = 20)
plot(u_rot[, 11:12], col = pop, pch = 20)
plot(u_rot[, 13:14], col = pop, pch = 20)
plot(u_rot[, 15:16], col = pop, pch = 20)
plot(u_rot[, 17:18], col = pop, pch = 20)
plot(u_rot[, 19:20], col = pop, pch = 20)
