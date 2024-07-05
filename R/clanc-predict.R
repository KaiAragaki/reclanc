#' Predict from a `clanc`
#'
#' @param object A `clanc` object.
#'
#' @param new_data A data frame or matrix of new predictors.
#'
#' @param type A single character. The type of predictions to generate.
#' Valid options are:
#'
#' - `"numeric"` for numeric predictions.
#'
#' @param ... Not used, but required for extensibility.
#'
#' @return
#'
#' A tibble of predictions. The number of rows in the tibble is guaranteed
#' to be the same as the number of rows in `new_data`.
#'
#' @examples
#' train <- mtcars[1:20,]
#' test <- mtcars[21:32, -1]
#'
#' # Fit
#' mod <- clanc(mpg ~ cyl + log(drat), train)
#'
#' # Predict, with preprocessing
#' predict(mod, test)
#'
#' @export
predict.clanc <- function(object, new_data, type, ...) {
  forged <- hardhat::forge(new_data, object$blueprint)
  rlang::arg_match(type, valid_clanc_predict_types())
  predict_clanc_bridge(type, object, forged$predictors)
}

valid_clanc_predict_types <- function() {
  c("numeric", "class")
}

# ------------------------------------------------------------------------------
# Bridge

predict_clanc_bridge <- function(type, model, predictors) {
  predictors <- t(predictors)

  predict_function <- get_clanc_predict_function(type)
  predictions <- predict_function(model, predictors)
  hardhat::validate_prediction_size(predictions, t(predictors))

  predictions
}

get_clanc_predict_function <- function(type) {
  switch(
    type,
    numeric = predict_clanc_numeric,
    class = predict_clanc_class
  )
}

# ------------------------------------------------------------------------------
# Implementation

predict_clanc_class <- function(model, predictors) {
  filter_predictors(model, predictors) |>
    calc_dist() |>
    predict_class() |>
    hardhat::spruce_class()
}

predict_clanc_numeric <- function(model, predictors) {

  predictions <- rep(1L, times = nrow(predictors))
  hardhat::spruce_numeric(predictions)
}
