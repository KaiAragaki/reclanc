#' @param expression a data.frame containing `sample_id` and `class`, at least
#' @param priors either a number, a numeric vector with length equal to the
#'   number of unique classes, "equal", or "class". If "equal", the same prior
#'   will be given to all classes. If "class", a prior based on class frequency
#'   will be used.
add_priors <- function(expression, priors) {
  priors <- validate_priors(expression, priors)

  df <- data.frame(
    class = levels(expression$class),
    n = as.numeric(table(expression$class))
  )

  if (is.numeric(priors)) df$prior <- priors
  if (priors == "equal") df$prior <- 1 / nrow(df)
  if (priors == "class") df$prior <- df$n / nrow(df)

  df <- dplyr::select(df, -.data$n)

  dplyr::left_join(expression, df, by = "class")
}

validate_priors <- function(expression, priors) {

  if (is.character(priors)) {
    stopifnot(length(priors) == 1)
    stopifnot(priors %in% c("class", "equal"))
  }

  if (is.numeric(priors)) {
    stopifnot(
      length(priors) == length(levels(expression$class)),
      sum(priors) == 1,
      all(priors >= 0)
    )
  }

}
