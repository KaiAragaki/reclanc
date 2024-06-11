pooled_sd_clanc <- function(X, ID) {
  m = nrow(X)
  n = ncol(X)
  p = ncol(ID)
  nn = drop(t(ID) %*% rep(1, n))

  mn = t(t(X %*% ID) / nn)
  dif2 = (X - mn %*% t(ID)) ^ 2

  sqrt(drop(dif2 %*% rep(1 / (n - p), n)))
}
