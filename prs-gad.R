library(bigsnpr)

plink <- download_plink()
celiac <- "../thesis-celiac/Dubois2010_data/FinnuncorrNLITUK3hap550"

system(glue::glue(
  "{plink} --score grs.txt --bfile {celiac}"
))


d <- read.table("plink.profile", header=TRUE)
intercept <- -0.757226
d$GRS <- d$SCORE * d$CNT + intercept
# write.table(d[, c("FID", "IID", "GRS")],
#             file="profile.txt", row.names=FALSE, col.names=FALSE)



ind.test500 <- ind.test[celiac$fam$family.ID[ind.test] != "UK1"]
preds500 <- preds2[celiac$fam$family.ID[ind.test] != "UK1"]

preds_gad <- d$GRS[match(celiac$fam$sample.ID[ind.test500], d$IID)]

qplot(preds500, preds_gad) +
  geom_smooth(method = "lm", color = "blue") +
  labs(x = "My pred (AUC=88.6%)", y = "Gad pred (AUC=89.3%)")

AUC(preds500, y[ind.test500])
AUC(preds_gad, y[ind.test500])
