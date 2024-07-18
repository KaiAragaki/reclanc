process_priors <- function(processed, active, priors, verbosity, ...) {
  stopifnot(is.character(priors) || is.numeric(priors))

  n_levels <- length(levels(processed$outcomes[[1]]))

  if (is.character(priors)) {
    stopifnot(length(priors) == 1, priors %in% c("class", "equal"))
    if (priors == "class") {
      freqs <- as.numeric(table(processed$outcomes))
      priors <- freqs / nrow(processed$outcomes)
    }
    if (priors == "equal") priors <- 1 / n_levels
  }

  if (length(priors) == 1) priors <- rep(priors, times = n_levels)
  if (length(priors) != n_levels) {
    cli::cli_abort(c(
      "Length of priors must be equal to 1 or number of unique class labels",
      "i" = "Length of priors: {length(priors)}",
      "i" = "Number of unique class labels: {n_levels}"
    ))
  }
  if (any(priors < 0)) {
    cli::cli_abort("All priors must be positive")
  }
  if (sum(priors) > 1) {
    cli::cli_abort("Sum of priors ({sum(priors)}) is greater than 1")
  }
  if (sum(priors) < 1 && verbosity %in% c("all", "warn")) {
    cli::cli_warn("Sum of priors ({sum(priors)}) is less than 1")
  }
  priors
}

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
