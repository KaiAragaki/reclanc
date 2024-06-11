#' @param data expression data. (m x n) matrix of class numeric.
#' @param id class IDs. n vector of class numeric.
#' @param priors either a numeric vector of length p, "equal", or "class". If
#'   "equal", equal probabilities will be used for each class. If "class", the
#'   proportions of each class in the training data will be used as the prior.
#' @param active how many active features to consider? can either be a single
#'   number or a vector containing a range of values to consider.
cvClanc <- function(data, id, priors = "equal", active = 1:10, folds = 5) {
  cvIdx <- make_balanced_folds(id, folds)
  nrow_data <- nrow(data)
  n <- ncol(data)
  p <- length(unique(id))
  nn <- as.numeric(table(id))
  folds <- length(cvIdx)

  class_priors <- make_class_priors(priors = priors, data = data, id = id)
  d <- length(active)

  cv.error <- array(rep(0, d * folds * p), dim = c(d, folds, p))
  cv.err.cnt.cls <- matrix(NA, nrow = d, ncol = p)
  cv.err.prpn.cls <- matrix(NA, nrow = d, ncol = p)
  n.features.cls <- matrix(NA, nrow = d, ncol = p)
  n.features.ttl <- rep(NA, d)

  ID <- model.matrix(~ factor(id) - 1) # one-hot encoded matrix
  dimnames(ID) <- list(NULL, names(nn))

  ## cross validation
  cat("CV:")
  for(i in seq_len(folds)) {
    cat(i)

    ## form initial statistics
    v <- length(cvIdx[[i]]) # n samples in this fold
    X <- data[, -cvIdx[[i]]] # samples not in this fold
    Y <- data[, cvIdx[[i]]] # samples in this fold
    jd <- id[-cvIdx[[i]]] # ids not in this fold
    truth <- id[cvIdx[[i]]] # ids in this fold
    JD <- ID[-cvIdx[[i]], ] #

    mm <- table(jd)
    m.k <- sqrt(1 / mm - 1 / (n - v))

    ## pooled standard deviations
    p.sd <- pooledSDClanc(X, JD)

    ## class- and overall-centroids
    cntrd.k <- scale(X %*% JD, center = FALSE, scale = mm)
    cntrd.o <- drop(X %*% rep(1 / (n - v), n - v))

    ## form statistics and order them
    d.k <- scale((cntrd.k - cntrd.o) / p.sd, center = FALSE, scale = m.k)
    d.k.ord <- orderStatsClanc(d.k = d.k)

    ## select genes, update inactive centroid components
    for(j in seq_len(d)) {
      # FIXME
      aa <- ifelse(is.matrix(active), active[j, ], active[j])
      selected <- selectClanc(d.k = d.k, d.k.ord = d.k.ord, active = aa)
      active.idx <- seq_len(nrow_data)[drop(selected %*% rep(1, p)) != 0]
      cntrds <- cntrd.k[active.idx, ]
      for(k in seq_len(p))
        cntrds[selected[active.idx, k] == 0, k] = cntrd.o[active.idx][selected[active.idx, k] == 0]

      ## classify test sample and assess error status
      for(k in seq_len(v)) {
        dd <- distClanc(data = Y[active.idx, k], cntrds = cntrds, sd = p.sd[active.idx], prior = class_priors)

        if(match(min(dd), dd) != truth[k])
          cv.error[j, i, truth[k]] = cv.error[j, i, truth[k]] + 1
      }
    }
  }

  ## record numbers and proportions of errors
  for(i in seq_len(p)) {
    if(d > 1)
      cv.err.cnt.cls[, i] <- apply(cv.error[, , i], 1, sum)
    else
      cv.err.cnt.cls[, i] <- sum(cv.error[, , i])
    cv.err.prpn.cls[, i] <- cv.err.cnt.cls[, i] / nn[i]
  }
  cv.err.cnt.ttl <- apply(cv.err.cnt.cls, 1, sum)
  cv.err.prpn.ttl <- cv.err.cnt.ttl / n

  cat("\n")
  list(classErrors = cv.err.prpn.cls, overallErrors = cv.err.prpn.ttl, prior = class_priors)
}
