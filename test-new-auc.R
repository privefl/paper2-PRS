# plot(c(0, 1), c(0, 1))
# bigstatsr::pasteLoc(8)

pts <- list(x = c(0, 0.007, 0.027, 0.034, 0.101, 0.196, 0.351, 0.601, 0.839, 1), 
            y = c(0, 0.134, 0.347, 0.545, 0.701, 0.789, 0.841, 0.895, 0.942, 1))
plot(pts, type = "b")
pts2 <- list(x = c(0, 0.02, 0.051, 0.073, 0.113, 0.152, 0.222, 0.404, 0.674, 1), 
             y = c(0, 0.097, 0.238, 0.369, 0.539, 0.698, 0.823, 0.916, 0.962, 1))
points(pts2, col = "red", type = "b")

f1 <- function(pts) {
  w <- diff(pts$x)
  h <- pts$y[-length(pts$y)] + diff(pts$y) / 2 
  drop(crossprod(h, w))
}
f1(pts)
f1(pts2)


f2 <- function(pts) {
  x <- pts$x
  y <- pts$y
  ind <- seq_len(length(x) - 1)
  parts <- sapply(ind, function(i) {
    d <- 1 - x[[i]]
    c <- 1 - x[[i + 1]]
    a <- y[[i]]
    b <- y[[i + 1]]
    (b - a) * (d^2 + c * d + c^2) / 3 + (a * d - b * c) * (d + c) / 2
  })
  sum(parts) * 2
}

f2(pts)
f2(pts2)

plot(1 - pts$x, pts$y)
curve(x + 0)
curve(1 / exp(10*(1-x)))
