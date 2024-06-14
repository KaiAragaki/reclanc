add_stats <- function(expression) {
  expression |>
    dplyr::mutate(mk = sqrt(1 / .data$class_size - 1 / .data$n_samples)) |>
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
