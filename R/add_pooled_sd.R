#' Calculate a pooled standard deviation for centroid
#'
#' @details This function will summarize the data, essentially removing
#'   sample-specific data
add_pooled_sd <- function(expression) {
  expression |>
    dplyr::mutate(
      n_samples = length(unique(.data$sample_id)),
      n_classes = length(levels(.data$class)),
      df = .data$n_samples - .data$n_classes
    ) |>
    dplyr::mutate(
      class_size = length(unique(.data$sample_id)),
      .by = "class"
    ) |>
    dplyr::mutate(
      sqd_error = (.data$expression - .data$class_centroid)^2,
      .by = c("class", "gene_id")
    ) |>
    dplyr::summarize(
      df = unique(.data$df),
      pooled_sd = sqrt(sum(.data$sqd_error, na.rm = TRUE) / .data$df),
      .by = c("class", "gene_id", "n_classes", "n_samples", "class_size",
              "class_centroid", "overall_centroid", "prior")
    )
}
