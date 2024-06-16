#' @param expression a matrix of expression data, where columns are genes and
#'   rows are samples.
#' @param classes factor vector of classes for each column in `data`
#' @param priors either a numeric vector of length p, "equal", or "class". If
#'   "equal", equal probabilities will be used for each class. If "class", the
#'   proportions of each class in the training data will be used as the prior.
#' @param active Number of genes to consider
cv_clanc <- function(expression,
                     classes,
                     priors = "equal",
                     active = 10,
                     n_folds = 5) {
  stopifnot(
    is.factor(classes),
    is.matrix(expression),
    ncol(expression) == length(classes)
  )
  expression <- tidy_expression(expression, classes) |>
    add_folds(n_folds) |>
    add_priors(priors)
  # cross validation
  out <- data.frame()
  for (i in seq_len(n_folds)) {
    test <- dplyr::filter(expression, .data$fold == i)
    train <-  dplyr::filter(expression, .data$fold != i)
    centroids <- make_centroids(train, active)
    new <- predict_from_centroids(test, centroids)
    out <- rbind(out, new)
  }
  out
}
