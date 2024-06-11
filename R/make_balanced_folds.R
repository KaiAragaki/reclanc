make_balanced_folds <- function(classes, nfolds) {
  stopifnot(nfolds <= max(table(classes)))
  if (nfolds > min(table(classes))) {
    warning(
      "More folds than unique items in class. ",
      "Some folds won't have all classes."
    )
  }

  n_classes <- length(unique(classes))
  class_idx <- split(seq(classes), classes)

  scrambled_folds <- class_idx |>
    # *Intra*-fold scramble. Size parameter needed in case class size = 1
    lapply(\(x) sample(x, size = length(x))) |>
    # Split as evenly as possible into nfolds
    lapply(
      \(x) split(x, cut(seq_along(x), nfolds, labels = FALSE))
    ) |>
    # *Inter*-fold scramble, so uneven numbers don't always go to one fold
    lapply(sample)

  # Extend each class to the nfolds length
  for (i in seq_along(scrambled_folds)) {
    diff <- nfolds - length(scrambled_folds[[i]])
    if (diff > 0) {
      scrambled_folds[[i]] <- c(
        scrambled_folds[[i]], rep(NA, diff)
      )
    }
  }

  # Append each fold of each class together
  res <- vector(mode = "list", length = nfolds)
  for (i in seq_len(nfolds)) {
    for (j in seq_len(n_classes)) {
      res[[i]] <- c(res[[i]], scrambled_folds[[j]][[i]])
    }
  }

  lapply(res, \(x) as.numeric(na.omit(x)))
}
