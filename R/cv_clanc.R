#' @param expression a matrix of expression data, where columns are genes and
#'   rows are samples.
#' @param classes factor vector of classes for each column in `data`
#' @param priors either a numeric vector of length p, "equal", or "class". If
#'   "equal", equal probabilities will be used for each class. If "class", the
#'   proportions of each class in the training data will be used as the prior.
#' @param active Number of genes to consider
cv_clanc <- function(expression,
                     classes,
                     priors = "equal",
                     active = 10,
                     n_folds = 5) {
  stopifnot(
    is.factor(classes),
    is.matrix(data),
    ncol(data) == length(classes)
  )
  expression <- tidy_expression(expression, classes) |>
    add_folds(n_folds) |>
    add_priors(priors)

  # cross validation
  for (i in seq_len(n_folds)) {
    current_fold <- dplyr::filter(expression, .data$fold == i)
    other_folds <-  dplyr::filter(expression, .data$fold != i)

    class_stats <- other_folds |>
      add_centroids() |>
      add_pooled_sd() |>
      add_stats()

    ## select genes, update inactive centroid components
    class_stats <- class_stats |>
      mark_active_genes(active) |>
      dplyr::filter(any(.data$gets), .by = "gene") |>
      dplyr::mutate(
        final_centroid = dplyr::if_else(
          .data$gets, .data$class_centroid, .data$overall_centroid
        )
      )

    current_fold |>
      dplyr::rename(truth = class) |>
      dplyr::semi_join(by = "gene_id") |>
      dplyr::full_join(class_stats, by = "gene_id") |>
      dplyr::mutate(
        dist = ((.data$expression - .data$class_centroid) / .data$pooled_sd)^2
      ) |>
      dplyr::summarize(
        sum_dist = sum(.data$dist) - 2 * log(.data$prior),
        .by = c("sample", "class", "truth")
      ) |>
      dplyr::filter(
        .data$sum_dist == min(.data$sum_dist),
        .by = "class"
      )
  }
}
