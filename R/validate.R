# Input validation
validate_summarized_experiment <- function() {

}


validate_se <- function(se, formula, priors, active) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  se <- validate_dv(se, formula, priors, active)
  se <- validate_iv(se, formula, priors, active)
  se
}

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

validate_formula <- function(formula) {
  stopifnot(inherits(formula, "formula"))
  if (length(all.vars(formula)) != 2) {
    cli::cli_abort(
      "formula must contain exactly one independent and one dependent variable"
    )
  }
  formula
}

validate_priors <- function(se, formula, priors) {
  stopifnot(is.character(priors) || is.numeric(priors))
  if (is.character(priors)) {
    stopifnot(
      length(priors) == 1,
      priors %in% c("class", "equal")
    )
  }

  dv <- get_var_info_from_form(se, fromula)
  n_levels <- length(levels(dv$data))

  if (length(priors) == 1) priors <- rep(priors, times = n_levels)

  if (is.numeric(priors)) {
    stopifnot(
      length(priors) == n_levels,
      sum(priors) == 1,
      all(priors >= 0)
    )
  }

  priors
}

validate_active <- function(se, formula, active) {
  stopifnot(is.numeric(active))
  dv <- get_var_info_from_form(se, formula)
  n_classes <- length(levels(dv$data))

  if (!length(active) %in% c(1, n_classes)) {
    cli::cli_abort(
      "{.code length(active)} must equal # levels in {.var {dv$name}} or 1",
      "i" = "{.code length(active)}: {.var {length(active)}}",
      "i" = "Levels in {.var {dv$name}}: {.var {n_classes}}"
    )
  }

  if (sum(active) > nrow(se)) {
    cli::cli_abort(
      "More active genes requested than exist in data",
      "i" = "Active genes requested for all classes: {.var {sum(active)}}",
      "i" = "Number in data: {.var {nrow(se)}}"
    )
  }

  active
}
