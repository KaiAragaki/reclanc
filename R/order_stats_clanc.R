#' @param dk a matrix of statistics, where each column is a class and each row
#'   is a gene
order_stats_clanc <- function(dk) {
  dk_ranks <- apply(dk, 2, \(x) order(abs(x), decreasing = TRUE))
  dk_ordered <- apply(dk, 2, \(x) sort(abs(x), decreasing = TRUE))
  dk_no_ties <- apply(dk_ordered, 2, \(x) !duplicated(x))
  list(d.k.rnks = dk_ranks, d.k.o = dk_ordered, no.ties = dk_no_ties)
}
