test_data <- function() {
  n_same <- 10000
  n_diff <- 100
  total <- n_same + n_diff
  p_same <- n_same / total
  p_diff <- n_diff / total
  probs <- c(p_same, p_diff)
  a_means <- sample(c(100, 110), total, probs, replace = TRUE)
  b_means <- sample(c(100, 100), total, probs, replace = TRUE)
  c_means <- sample(c(100, 200), total, probs, replace = TRUE)
  d_means <- sample(c(100, 150), total, probs, replace = TRUE)
  e_means <- sample(c(100, 10), total, c(p_same, 20 / 10000), replace = TRUE)
  n_classes <- rep(1:5, each = 20)
  mean_matrix <- matrix(
    c(a_means, b_means, c_means, d_means, e_means),
    nrow = total
  )
  apply(mean_matrix[, n_classes], c(1, 2), \(x) rpois(1, lambda = x))
}

test_dk <- function() {
  test_dk <- expand.grid(
    gene_id = 1:10,
    class = factor(letters[1:3])
  )
  test_dk$stat <- c(-10, 7, 5, 5, 3, 1.5, 2, 23, 8, 21,
                    1:10,
                    100, 99, 98, 97, -101, 3, 33, 333, 2, 109)
  test_dk |>
    dplyr::mutate(
      abs_stat = abs(.data$stat),
      rank = rank(-.data$abs_stat, ties.method = "min"), .by = class
    )
}

test_active <- function() {
  data.frame(
    class = factor(letters[1:3]),
    active = c(4, 2, 3)
  )
}
