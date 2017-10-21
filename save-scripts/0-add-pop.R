library(bigsnpr)

celiac <- snp_attach("backingfiles/celiacQC.rds")
G <- celiac$genotypes  # 15155 x 281122 (would take 32 Gb)

# get population
pop.files <- list.files(path = "../thesis-celiac/Dubois2010_data/", 
                        pattern = "cluster_*", full.names = TRUE)
pop <- snp_getSampleInfos(celiac, pop.files)[[1]]
pop.names <- c("Netherlands", "Italy", "UK1", "UK2", "Finland")
celiac$fam$family.ID <- pop.names[pop]
celiac <- snp_save(celiac)
celiac$fam$family.ID
