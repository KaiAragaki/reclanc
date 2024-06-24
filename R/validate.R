validate_dv <- function(se, formula, priors, active) {
  dv <- get_var_info_from_form(se, formula, 1)

  if (is.null(dv$data)) {
    cli::cli_abort(
      "{.var {dv$name}} does not exist in {.code colData(se)}"
    )
  }

  single_prior_and_active <-
    (length(unique(priors)) == 1) && (length(unique(active)) == 1)

  if (!is.factor(dv$data) && single_prior_and_active) {
    cli::cli_inform("{.var {dv$name}} is not a factor, converting")
    se[[dv$name]] <- factor(dv$data)
  }

  if (!is.factor(dv$data) && !single_prior_and_active) {
    cli::cli_abort(c(
      "{.var {dv$name}} is not a factor and length(prior/active) > 1",
      "i" = "Mapping from class to prior/active is ambiguous",
      "i" = "Convert {.var {dv$name}} to factor and try again."
    ))
  }

  se
}

validate_iv <- function(se, formula) {
  iv <- get_var_info_from_form(se, formula, 2)
  if (!iv$name %in% colnames(rowData(se))) {
    cli::cli_abort(
      "{.var {iv$name}} does not exist in {.code rowData(se)}"
    )
  }

  if (any(is.na(iv$data))) {
    cli::cli_warn("{.var {iv_name}} contains NA values. Removing them.")
    se <- se[which(is.na(iv$data)), ]
  }

  if (any(duplicated(iv$data))) {
    cli::cli_abort("Duplicate values found in {.var iv_name}")
  }

  rownames(se) <- iv$data

  se
}

process_priors <- function(processed, active, priors, ...) {
  stopifnot(is.character(priors) || is.numeric(priors))

  n_levels <- length(levels(processed$outcomes$.outcome))

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
  if (sum(priors) < 1) {
    cli::cli_warn("Sum of priors ({sum(priors)}) is less than 1")
  }
  priors
}

process_active <- function(processed, active, priors, ...) {
  stopifnot(is.numeric(active))

  n_levels <- length(levels(processed$outcomes$.outcome))

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

validate_factors <- function() {

}
