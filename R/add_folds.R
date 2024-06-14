add_folds <- function(expression, n_folds) {
  validate_folds(expression, n_folds)

  fold_table <- make_fold_table(expression, n_folds)
  dplyr::left_join(expression, fold_table, by = "sample_id")
}

validate_folds <- function(expression, n_folds) {
  stopifnot(is.factor(expression$class))

  if (n_folds > min(table(expression$class))) {
    stop("More folds than unique items in class")
  }
}

#' @param expression tidy expression, usually output from `tidy_expression`
#' @param n_folds number of desired folds
#' @return A data.frame with columns `class`, `fold`, and `cut`. The table
#'   includes all possible class and fold combinations, with `cut` randomly
#'   assigned to each. See details.
#' @details The classes are divided as evenly as possible into n_folds 'cuts'.
#'   As these cuts are not necessarily even (dividing 10 samples 3 ways, for
#'   instance), which cut goes to which fold is scrambled so that one fold isn't
#'   always getting the 'large' or 'small' cuts
make_fold_table <- function(expression, n_folds) {
  fold_table <- expression |>
    dplyr::select(.data$sample_id, .data$class) |>
    dplyr::distinct() |>
    # Shuffle all samples within class in case there's an inherent order
    dplyr::slice_sample(prop = 1, by = "class") |>
    dplyr::mutate(
      cut = as.numeric(cut(seq_len(dplyr::n()), n_folds)),
      .by = "class"
    )

  classes <- levels(expression$class)
  cut_class_fold_map <- expand.grid(
    cut = seq_len(n_folds),
    class = classes
  ) |>
    dplyr::mutate(
      fold = replicate(
        length(classes), sample(n_folds), simplify = FALSE
      ) |>
        unlist()
    )

  dplyr::left_join(fold_table, cut_class_fold_map, by = c("class", "cut")) |>
    dplyr::select(.data$sample_id, .data$fold)
}
