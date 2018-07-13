add_sens_FDP2 <- function(x, corr) {
  
  bigstatsr:::assert_class(corr, "dsCMatrix")
  corr <- as(abs(corr), "dgTMatrix")
  
  alt <- map2(x$true_set, x$set, ~ {
    if (length(.y) == 0) return(c(0, 0))
    corr_set <- corr[.x, .y, drop = FALSE]
    AltSens <-     sum(tapply(corr_set@x, corr_set@i, max)) / nrow(corr_set)
    AltFDP  <- 1 - sum(tapply(corr_set@x, corr_set@j, max)) / ncol(corr_set)
    c(AltSens, AltFDP)
  }) %>%
    transpose()
  
  x %>%
    mutate(AltSens = unlist(alt[[1]]),
           AltFDP  = unlist(alt[[2]]))
}

add_sens_FDP3 <- function(x, corr) {
  
  bigstatsr:::assert_class(corr, "dsCMatrix")
  corr <- as(abs(corr), "dgTMatrix")
  
  alt <- map2(x$true_set, x$set, ~ {
    if (length(.y) == 0) return(c(0, 0))
    bool.x <- bool.y <- logical(nrow(corr))
    bool.x[.x] <- TRUE
    bool.y[.y] <- TRUE
    ind <- which_cond(corr@i, corr@j, bool.x, bool.y)
    AltSens <-     sum(tapply(corr@x[ind], corr@i[ind], max)) / length(.x)
    AltFDP  <- 1 - sum(tapply(corr@x[ind], corr@j[ind], max)) / length(.y)
    c(AltSens, AltFDP)
  }) %>%
    transpose()
  
  x %>%
    mutate(AltSens = unlist(alt[[1]]),
           AltFDP  = unlist(alt[[2]]))
}

x <- list.files("results5", full.names = TRUE)[1:20] %>%
  map_dfr(~readRDS(.x)) %>%
  select(true_set, set) %>%
  as_tibble()
  
corr <- readRDS("backingfiles/corr2.rds")

# system.time(print(add_sens_FDP(x, corr)))
# system.time(print(add_sens_FDP2(x, corr)))
Rcpp::sourceCpp('save-scripts/compute-alt-metrics.cpp')
system.time(print(add_sens_FDP3(x, corr)))

library(Matrix)
X <- abs(corr[1:100, 1:100 + 100])
system.time(
  tmp <- apply(X, 2, max)
)
sum(tmp)

X2 <- as(X, "dgTMatrix")
system.time({
  X3 <- data_frame(i = X2@i, j = X2@j, x = X2@x)
  tmp2 <- group_by(X3, j) %>%
    summarise(max(x))
})
sum(tmp2$`max(x)`)

system.time(
  tmp3 <- tapply(X2@x, X2@j, max)
)
sum(tmp3)

rbind(tmp[as.integer(names(tmp3)) + 1], tmp3)
