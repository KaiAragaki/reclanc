process_priors <- function(processed, active, priors) {
  stopifnot(is.character(priors) || is.numeric(priors))
  classes <- processed$outcomes[[1]]
  priors <- process_priors_length(priors, length(levels(classes)))
  if (is.character(priors)) return(process_priors_char(priors, classes))
  if (is.numeric(priors)) return(process_priors_num(priors))
}

process_priors_length <- function(priors, n_levels) {
  if (is.character(priors)) {
    if (length(priors) != 1) {
      cli::cli_abort("Length of priors must be 1 if priors is a character")
    }
  }
  if (length(priors) == 1) priors <- rep(priors, times = n_levels)
  if (length(priors) == n_levels) {
    cli::cli_abort(c(
      "Length of priors must be equal to 1 or number of unique class labels",
      "i" = "Length of priors: {length(priors)}",
      "i" = "Number of unique class labels: {n_levels}"
    ))
  }
}

# Assumes input is correct length
# eg a length of 1 has been recycled
process_priors_num <- function(priors) {
  if (any(priors < 0)) cli::cli_abort("All priors must be positive")
  # HACK Should probably use machine epsilon
  if (sum(priors) > 1.01) cli::cli_abort("Sum of priors > 1 ({sum(priors)})")
  if (sum(priors) < 0.99)
    cli::cli_warn("Sum of priors < 1 ({sum(priors)})")
  priors
}

# Assumes priors input is length 1
process_priors_char <- function(priors, classes) {
  if (priors == "class") {
    freqs <- as.numeric(table(classes))
    return(freqs / length(classes))
  }
  if (priors == "equal") {
    n_class <- length(levels(classes))
    return(rep(1 / n_class, n_class))
  }
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
