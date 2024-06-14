add_stats <- function(expression) {
  n_samples <- length(unique(expression$sample_id))
  expression |>
    dplyr::mutate(mk = sqrt(1 / dplyr::n() - 1 / n_samples), .by = "class") |>
    dplyr::mutate(
      stat = ((.data$class_centroid - .data$overall_centroid) / .data$pooled_sd) * .data$mk,
      abs_stat = abs(.data$stat)
    ) |>
    dplyr::mutate(
      class_rank = rank(-.data$abs_stat, ties.method = "min"),
      no_tie = !duplicated(.data$class_rank),
      .by = "class"
    ) |>
    dplyr::select(-"mk")
}
