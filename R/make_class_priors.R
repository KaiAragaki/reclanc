make_class_priors <- function(priors, id, data) {
  n_classes <- length(unique(id))
  n_samples <- ncol(data)
  samples_per_class <- as.numeric(table(id))
  class_proportions <- n_samples / samples_per_class

  if (is.numeric(priors)) {
    stopifnot(length(priors) == n_classes, sum(priors) == 1, all(priors >= 0))
    return(priors)
  }

  stopifnot(priors %in% c("equal", "class"))
  if (priors == "equal") return(rep(1 / n_classes, n_classes))
  if (priors == "class") return(class_proportions)
}
