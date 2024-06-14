#' @param expression a matrix of expression, where columns are samples and rows
#'   are genes.
#' @param classes a factor vector of classes for each column in `expression`
#' @return A data.frame with columns `expression`, `class`, `sample_id`, and
#'   `gene_id`
tidy_expression <- function(expression, classes) {
  stopifnot(is.matrix(expression), is.factor(classes))

  if (is.null(rownames(expression))) {
    rownames(expression) <- seq_len(nrow(expression))
  }

  if (is.null(colnames(expression))) {
    colnames(expression) <- seq_len(ncol(expression))
  }

  out <- expression |>
    matrix(ncol = 1) |>
    as.data.frame() |>
    cbind(rep(classes, each = nrow(expression))) |>
    cbind(rep(colnames(expression), each = nrow(expression))) |>
    cbind(rownames(expression))
  colnames(out) <- c("expression", "class", "sample_id", "gene_id")
  out
}
