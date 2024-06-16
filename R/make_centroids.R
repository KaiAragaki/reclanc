#' @param expression tidy expression - where each row represents a single
#'   expression for a single sample for a single gene.
#' @param active the number of genes to select for each centroid (can be a
#'   vector to specify differing numbers of genes per class)
make_centroids <- function(expression, active) {
  expression |>
    add_class_and_overall_means() |>
    add_pooled_sd() |>
    add_stats() |>
    mark_active_genes(active) |>
    dplyr::filter(any(.data$gets), .by = "gene_id") |>
    dplyr::mutate(
      final_centroid = dplyr::if_else(
        .data$gets, .data$class_centroid, .data$overall_centroid
      )
    )
}
