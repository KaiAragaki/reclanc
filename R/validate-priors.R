process_priors <- function(processed, priors) {
  stopifnot(is.character(priors) || is.numeric(priors))
  classes <- processed$outcomes[[1]]
  if (is.character(priors)) return(process_priors_char(priors, classes))
  if (is.numeric(priors)) return(process_priors_num(priors, classes))
}

# Assumes input is correct length
# eg a length of 1 has been recycled
process_priors_num <- function(priors, classes) {
  priors <- process_priors_length(priors, classes)
  if (any(priors < 0)) cli::cli_abort("All priors must be positive")
  # HACK Should probably use machine epsilon
  if (sum(priors) > 1.01) cli::cli_abort("Sum of priors > 1 ({sum(priors)})")
  if (sum(priors) < 0.99)
    cli::cli_warn("Sum of priors < 1 ({sum(priors)})")
  priors
}

# Assumes priors input is length 1
process_priors_char <- function(priors, classes) {
  if (length(priors) != 1)
    cli::cli_abort("Length of priors must be 1 if priors is a character")

  if (priors == "class") {
    freqs <- as.numeric(table(classes))
    return(freqs / length(classes))
  }
  if (priors == "equal") {
    n_class <- length(levels(classes))
    return(rep(1 / n_class, n_class))
  }
}

process_priors_length <- function(priors, classes) {
  n_levels <- length(levels(classes))
  if (length(priors) == 1) priors <- rep(priors, times = n_levels)
  if (length(priors) != n_levels) {
    cli::cli_abort(c(
      "Length of priors must be equal to 1 or number of unique class labels",
      "i" = "Length of priors: {length(priors)}",
      "i" = "Number of unique class labels: {n_levels}"
    ))
  }
  priors
}
