#' @param data expression data. (m x n) matrix of class numeric.
#' @param id class IDs. n vector of class numeric.
#' @param priors either a numeric vector of length p, "equal", or "class". If
#'   "equal", equal probabilities will be used for each class. If "class", the
#'   proportions of each class in the training data will be used as the prior.
#' @param active how many active features to consider? can either be a single
#'   number or a vector containing a range of values to consider.
cv_clanc <- function(data, id, priors = "equal", active = 1:10, n_folds = 5) {
  fold_indices <- make_balanced_n_folds(id, n_folds)
  nrow_data <- nrow(data)
  n <- ncol(data)
  ncol_data <- ncol(data)
  p <- length(unique(id))
  n_classes <- length(unique(id))
  nn <- as.numeric(table(id))
  n_folds <- length(fold_indices)

  class_priors <- make_class_priors(priors = priors, data = data, id = id)
  d <- length(active)

  # Initializations
  cv.error <- array(rep(0, d * n_folds * n_classes), dim = c(d, n_folds, n_classes))
  cv.err.cnt.cls <- matrix(NA, nrow = d, ncol = n_classes)
  cv.err.prpn.cls <- matrix(NA, nrow = d, ncol = n_classes)
  n.features.cls <- matrix(NA, nrow = d, ncol = n_classes)
  n.features.ttl <- rep(NA, d)
  ID <- model.matrix(~ factor(id) - 1) # one-hot encoded matrix
  dimnames(ID) <- list(NULL, names(nn))

  ## cross validation
  cat("CV:")
  for(i in seq_len(n_folds)) {
    cat(i)

    ## form initial statistics
    v <- length(fold_indices[[i]]) # n samples in this fold
    X <- data[, -fold_indices[[i]]] # samples not in this fold
    Y <- data[, fold_indices[[i]]] # samples in this fold
    jd <- id[-fold_indices[[i]]] # ids not in this fold
    truth <- id[fold_indices[[i]]] # ids in this fold
    JD <- ID[-fold_indices[[i]], ] # Samples not in this fold (one-hot)

    mm <- table(jd)
    m.k <- sqrt(1 / mm - 1 / (ncol_data - v))

    ## pooled standard deviations
    p.sd <- pooled_sd_clanc(X, JD)

    ## class- and overall-centroids
    cntrd.k <- scale(X %*% JD, center = FALSE, scale = mm)
    cntrd.o <- drop(X %*% rep(1 / (ncol_data - v), ncol_data - v))

    ## form statistics and order them
    d.k <- scale((cntrd.k - cntrd.o) / p.sd, center = FALSE, scale = m.k)
    d.k.ord <- order_stats_clanc(d.k = d.k)

    ## select genes, update inactive centroid components

    # For each item in the active vector
    for (j in seq_len(d)) {
      aa <- active[j] # Number of active genes to consider
      selected <- select_clanc(d.k = d.k, d.k.ord = d.k.ord, active = aa)
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
