ttrees <- function(TTree, file.base, pheno.all, ind.train, ind.test,
                   n.trees = 100) {
  
  file0.jdb  <- paste0(file.base, ".jdb")
  file0.bloc <- paste0(file.base, ".bloc")
  
  # Write jdb file with new pheno
  tmpfile <- tempfile()
  write(c(rep("", 5), pheno.all), ncolumns = 1, tmpfile)
  # https://stackoverflow.com/a/7846550/6103040
  system(sprintf("awk 'FNR==NR{a[NR]=$1;next}{$2=a[FNR]}1' %s %s > %s",
                 tmpfile, file0.jdb, file.jdb <- paste0(tmpfile, ".jdb")))
  # Write indices of learning and validation sets
  file.learn <- paste0(tmpfile, "_learn.txt")
  file.val   <- paste0(tmpfile, "_val.txt")
  cat(ind.train - 1, file = file.learn, sep = "\t")
  cat(ind.test - 1,  file = file.val,   sep = "\t")
  
  timing <- system.time(
    system(glue::glue(
      "{TTree}",
      " -j {file.jdb}",
      " -m 3",
      " -b {file0.bloc}",
      " -l {file.learn}",
      " -v {file.val}", 
      " -t {n.trees} -k 1000 -c 5 -n 2000",
      " -x -s"
    ))
  )[3]
  
  file.roc <- sprintf("%s_k1000_m3_t%d_ic5_nmin2000_0000.roc",
                      file.jdb, n.trees)
  file.vim <- sub("\\.roc$", ".vim", file.roc)
  preds <- read.table(file.roc, header = FALSE)
  
  tibble(
    method   = "T-Trees",
    eval     = list(preds[match(ind.test - 1, preds[[1]]), 2:3]),
    timing   = timing,
    nb.preds = sum(read.table(file.vim)$V2 != 0)
  )
}

TTree <- "../../TTree-source/TTree"

test3 <- ttrees(TTree, "ttrees", pheno.all, ind.train, ind.test, n.trees = 100)
test4 <- ttrees(TTree, "ttrees", pheno.all, ind.train, ind.test, n.trees = 1e3)

str(test3)
str(test4)
plot(test3$eval$V2, test4$eval$V2)
AUC(test3$eval$V2, test3$eval$V3)
AUC(test4$eval$V2, test4$eval$V3)
