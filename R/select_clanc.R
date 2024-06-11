#' @param d.kiu
#' @param d.k.ord
#' @param active
select_clanc <- function(d.k, d.k.ord, active) {
  m = nrow(d.k)
  p = ncol(d.k)
  d.k.rnks = d.k.ord$d.k.rnks
  d.k.o = d.k.ord$d.k.o
  no.ties = d.k.ord$no.ties

  if(length(active) != p)
      active = rep(active, p)
  if(any(active >= m))
      selected = matrix(1, nrow = m, ncol = p)
  else {
    if(any((avail <- t(no.ties) %*% rep(1, m)) < active)) {
      active[avail < active] = avail[avail < active]

      warning("Not enough unique genes in at least one class.")
    }

    selected = matrix(0, nrow = m, ncol = p)

    jdx = vector("list", p)
    delta = rep(NA, p)
    for(i in 1:p) {
      jdx[[i]] = d.k.rnks[no.ties[, i], i][1:active[i]]
      delta[i] = d.k.o[no.ties[, i], i][active[i] + 1]
    }
    cls = rep(1:p, active)

    while(any((kdx <- table(unlist(jdx))) > 1)) {
      nms = as.numeric(names(kdx))
      for(j in as.numeric(names(kdx))[kdx > 1]) {
        ldx = cls[unlist(jdx) == j]

        ii = match(max(abs(d.k[j, ldx])), abs(d.k[j, ldx]))

        for(k in ldx[-ii]) {
          active[k] = active[k] + 1
          jdx[[k]] = c(jdx[[k]][-match(j, jdx[[k]])], d.k.rnks[no.ties[, k], k][active[k]])
          delta[k] = d.k.o[no.ties[, k], k][active[k] + 1]
        }
      }
    }

    for(i in 1:p)
      selected[jdx[[i]], i] = 1
  }

  return(selected)
}
