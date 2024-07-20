process_active <- function(processed, active, priors, ...) {
  stopifnot(is.numeric(active))

  n_levels <- length(levels(processed$outcomes[[1]]))

  if (length(active) == 1) active <- rep(active, n_levels)

  if (length(active) != n_levels) {
    cli::cli_abort(c(
      "Length of `active` must be equal to 1 or number of unique class labels",
      "i" = "Length of `active`: {length(active)}",
      "i" = "Number of unique class labels: {n_levels}"
    ))
  }

  if (sum(active) > ncol(processed$predictors)) {
    cli::cli_abort(
      "More active genes requested than exist in data",
      "i" = "Active genes requested for all classes: {sum(active)}",
      "i" = "Number in data: {ncol(processed$predictors)}"
    )
  }

  active
}

validate_classes <- function(processed, active, priors, ...) {
  classes <- processed$outcomes[[1]]

  if (!is.factor(classes)) cli::cli_abort("Classes must be a factor")

  n_classes <- length(levels(classes))

  if (n_classes == 1)
    cli::cli_abort("Finding centroids with only one class doesn't make sense.")
}
