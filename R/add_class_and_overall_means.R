add_class_and_overall_means <- function(expression) {
  expression |>
    dplyr::mutate(
      class_centroid = mean(.data$expression, na.rm = TRUE),
      .by = c("class", "gene_id")
    ) |>
    dplyr::mutate(
      overall_centroid = mean(.data$expression, na.rm = TRUE),
      .by = "gene_id"
    )
}
