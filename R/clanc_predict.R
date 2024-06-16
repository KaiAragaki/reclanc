clanc_predict <- function(expression, classes, centroids) {
  expression <- tidy_expression(expression, classes)
  # TODO Check to see if expression has active genes
  predict_from_centroids(expression, centroids)
}

# Internal version
predict_from_centroids <- function(expression, centroids) {
  expression |>
    dplyr::rename(truth = .data$class) |>
    dplyr::select(-"prior") |>
    dplyr::semi_join(centroids, by = "gene_id") |>
    dplyr::full_join(centroids, by = "gene_id") |>
    dplyr::mutate(
      dist = ((.data$expression - .data$class_centroid) / .data$pooled_sd)^2
    ) |>
    dplyr::summarize(
      sum_dist = sum(.data$dist),
      .by = c("sample_id", "class", "truth", "prior")
    ) |>
    dplyr::mutate(final_dist = .data$sum_dist - 2 * log(.data$prior)) |>
    dplyr::filter(
      .data$final_dist == min(.data$final_dist),
      .by = "sample_id"
    )
}
