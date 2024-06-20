# This function was ripped from {recipes}:
# https://github.com/tidymodels/recipes

# Some assumptions have been made to reduce the amount needed to grok
# Assumes left and right side size are always greater than 0
# Assumes input is always formula (will be enforced upstream)
get_rhs_vars <- function(formula, data, no_lhs = FALSE) {

  if (no_lhs) {
    formula <- rlang::new_formula(lhs = NULL, rhs = rlang::f_rhs(formula))
  }

  outcomes_names <- all.names(
    rlang::f_lhs(formula),
    functions = FALSE,
    unique = TRUE
  )

  predictors_names <- all.names(
    rlang::f_rhs(formula),
    functions = FALSE,
    unique = TRUE
  )

  if (any(predictors_names == ".")) {
    predictors_names <- predictors_names[predictors_names != "."]
    predictors_names <- c(predictors_names, colnames(data))
    predictors_names <- unique(predictors_names)
  }

  setdiff(predictors_names, outcomes_names)
}

form2args <- function(formula, data, ..., call = rlang::caller_env()) {
  ## check for in-line formulas
  inline_check(formula, data, call)

  ## use rlang to get both sides of the formula
  outcomes <- get_lhs_vars(formula, data)
  predictors <- get_rhs_vars(formula, data, no_lhs = TRUE)

  ## if . was used on the rhs, subtract out the outcomes
  predictors <- predictors[!(predictors %in% outcomes)]

  ## get `vars` from lhs and rhs of formula
  vars <- c(predictors, outcomes)

  ## subset data columns
  data <- data[, vars]

  list(data = data, predictors = predictors, outcomes = outcomes)
}

inline_check <- function(x, data, call) {
  funs <- fun_calls(x, data)
  funs <- funs[!(funs %in% c("~", "+", "-", "."))]

  if (length(funs) > 0) {
    cli::cli_abort(c(
      x = "Misspelled variable name or in-line functions detected.",
      i = "{cli::qty(length(funs))}The following function{?s}/misspelling{?s} \\
          {?was/were} found: {.and {.code {funs}}}.",
      i = "Use steps to do transformations instead.",
      i = "If your modeling engine uses special terms in formulas, pass \\
          that formula to workflows as a \\
          {.help [model formula](parsnip::model_formula)}."
    ), call = call)
  }

  invisible(x)
}

fun_calls <- function(f, data) {
  setdiff(all.names(f), colnames(data))
}

too_many_case_weights <- function(x, call = rlang::caller_env()) {
  n <- length(x)

  cli::cli_abort(
    c(
      "!" = "There should only be a single column with the role \\
      {.code case_weights}.",
      "i" = "In these data, there are {n} columns: {.var {x}}."
    ),
    call = call
  )
}

get_lhs_vars <- function(formula, data) {
  if (!rlang::is_formula(formula)) {
    formula <- as.formula(formula)
  }
  ## Want to make sure that multiple outcomes can be expressed as
  ## additions with no cbind business and that `.` works too (maybe)
  new_formula <- rlang::new_formula(lhs = NULL, rhs = rlang::f_lhs(formula))
  get_rhs_vars(new_formula, data)
}
