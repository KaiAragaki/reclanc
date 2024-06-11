#' @param exp An expression matrix, where rows are genes and columns are
#'   samples
#' @param classes_one_hot A one-hot encoded matrix of classes
#' @return A vector containing a pooled SD for each gene
pooled_sd_clanc <- function(exp, classes_one_hot) {
  df <- ncol(exp) - ncol(classes_one_hot)
  samples_per_class <- colSums(classes_one_hot)
  avg_gene_exp_per_class <- t(t(exp %*% classes_one_hot) / samples_per_class)
  squared_error <- (exp - (avg_gene_exp_per_class %*% t(classes_one_hot)))^2
  sqrt(rowSums(squared_error / df))
}
