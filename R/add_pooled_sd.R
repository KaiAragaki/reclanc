add_pooled_sd <- function(expression) {
  n_classes <- length(levels(expression$classes))
  n_samples <- length(unique(expression$sample_id))
  df <- n_samples - n_classes
  expression |>
    dplyr::mutate(
      class_mean_exp = mean(.data$expression, na.rm = TRUE),
      sqd_error = (.data$expression - .data$class_mean_exp)^2,
      .by = c("class", "gene_id")
    ) |>
    dplyr::mutate(
      pooled_sd = sqrt(sum(.data$sqd_error, na.rm = TRUE) / df),
      .by = "gene_id"
    ) |>
    dplyr::select(-c("class_mean_exp", "squared_error"))
}
