## Don't forget to change the phenotype
celiac2$fam$affection <- pheno
tmpfile <- tempfile()
snp_writeBed(celiac2, bedfile = paste0(tmpfile, ".bed"), ind.row = ind.train)


crossval <- "../../SparSNP/crossval.sh"
Sys.setenv(PATH = paste(normalizePath("../../SparSNP"), 
                        Sys.getenv("PATH"),
                        sep = ":"))

unlink("discovery/", recursive = TRUE, force = TRUE)
system(glue::glue("crossval.sh {tmpfile} sqrhinge 1"))

