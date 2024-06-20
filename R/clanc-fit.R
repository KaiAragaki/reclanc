#' Fit a `clanc`
#'
#' `clanc()` fits a model.
#'
#' @param x Depending on the context:
#'
#'   * A __data frame__ of predictors.
#'   * A __matrix__ of predictors.
#'   * A __recipe__ specifying a set of preprocessing steps
#'     created from [recipes::recipe()].
#'
#' @param y When `x` is a __data frame__ or __matrix__, `y` is the outcome
#' specified as:
#'
#'   * A __data frame__ with 1 numeric column.
#'   * A __matrix__ with 1 numeric column.
#'   * A numeric __vector__.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both the predictors and the outcome.
#'
#' @param formula A formula specifying the outcome terms on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @return
#'
#' A `clanc` object.
#'
#' @examples
#' predictors <- mtcars[, -1]
#' outcome <- mtcars[, 1]
#'
#' # XY interface
#' mod <- clanc(predictors, outcome)
#'
#' # Formula interface
#' mod2 <- clanc(mpg ~ ., mtcars)
#'
#' # Recipes interface
#' library(recipes)
#' rec <- recipe(mpg ~ ., mtcars)
#' rec <- step_log(rec, disp)
#' mod3 <- clanc(rec, mtcars)
#'
#' @export
clanc <- function(x, ...) {
  UseMethod("clanc")
}

#' @export
#' @rdname clanc
clanc.default <- function(x, ...) {
  stop("`clanc()` is not defined for a '", class(x)[1], "'.", call. = FALSE)
}

# XY method - data frame

#' @export
#' @rdname clanc
clanc.data.frame <- function(x, y, active, priors, ...) {
  processed <- hardhat::mold(x, y)
  clanc_bridge(processed, active, priors, ...)
}

# XY method - matrix

#' @export
#' @rdname clanc
clanc.matrix <- function(x, y, active, priors, ...) {
  processed <- hardhat::mold(x, y)
  clanc_bridge(processed, active, priors, ...)
}

# Formula method

#' @export
#' @rdname clanc
clanc.formula <- function(formula, data, active, priors, ...) {
  # Oftentimes the formula will be incredibly wide.
  # R hates this.
  #
  # We get around this using some recent functions from recipes and pretend a
  # data.frame was supplied
  args <- form2args(formula, data)
  data <- args$data
  pred <- data[, which(colnames(data) %in% args$predictors)]
  outcomes <- data[, which(colnames(data) %in% args$outcomes)]

  clanc.data.frame(pred, outcomes, active, priors, ...)
}

# Recipe method

#' @export
#' @rdname clanc
clanc.recipe <- function(x, data, active, priors, ...) {
  processed <- hardhat::mold(x, data)
  clanc_bridge(processed, ...)
}

# ------------------------------------------------------------------------------
# Bridge

clanc_bridge <- function(processed, active, priors, ...) {
  expression <- processed$predictors
  classes <- processed$outcomes[[1]]
  class_data <- data.frame(class = sort(unique(classes)), active = active, priors = priors)
  fit <- clanc_impl(expression, class_data, classes)

  new_clanc(
    centroids = fit$centroids,
    blueprint = processed$blueprint
  )
}


# ------------------------------------------------------------------------------
# Implementation


# Expression - ***sample*gene*** matrix with names for each
#
# Classes + Active + Prior - data.frame
#
# Class to sample map - data.frame.
# colnames(class_sample_map) == c("sample", "class")
# Assumed to be arranged in the same order as expression rows
#
clanc_impl <- function(expression, class_data, classes) {
  fit <- clanc_fit(expression, class_data, classes)
  list(centroids = fit$centroids)
}

# TODO Allow things like "class" and "equal" to be provided to prior
