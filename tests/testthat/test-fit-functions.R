test_that("making ptypes works", {
  fit <- list(centroids = data.frame(gene = letters))
  processed <- list(
    blueprint = list(
      ptypes = list(
        outcomes = data.frame(.outcome = factor())
      )
    )
  )
  new_ptypes <- make_new_ptypes(fit, processed)
  expect_equal(colnames(new_ptypes$predictors), letters)
  expect_true(all(sapply(new_ptypes$predictors, class) ==  "numeric"))
  expect_true(all(sapply(new_ptypes, nrow) == 0))
})

test_that("class pooled sd calculation works", {
  n_class <- 4
  n_gene <- 3
  n_rep <- 3
  classes <- factor(rep(letters[seq_len(n_class)], each = n_rep))
  exp <- matrix(seq_len(n_class * n_gene * n_rep), ncol = n_gene)
  class_data <- data.frame(class = letters[seq_len(n_class)])
  expect_equal(
    calculate_pooled_sd(exp, class_data, classes),
    c(1, 1, 1)
  )
})


# TODO:
# Case of tie
# Case of tie and one active is full
# Case of more active than exist? Or is that handled upstream?
test_that("gene selection works", {

  abs_stats <- matrix(
    c(1, 2, 3,
      2, 3, 1,
      3, 1, 2),
    nrow = 3, byrow = TRUE,
    dimnames = list(letters[4:6], letters[1:3])
  )

  class_data <- data.frame(
    class = factor(letters[4:6]),
    active = rep(1, 3)
  )

  ranks <- t(apply(abs_stats, 1, \(x) rank(-x, ties.method = "min")))
  ties <- t(apply(ranks, 1, duplicated))
  selected <- select_genes(abs_stats, ranks, ties, class_data)
  expect_equal(selected, data.frame(class = factor(letters[4:6]), gene = letters[3:1]))
})


# TODO: Make sure this function respects factor order
test_that("class statistics calculation works", {
  n_class <- 4
  n_gene <- 3
  n_rep <- 3
  classes <- factor(rep(letters[seq_len(n_class)], each = n_rep))
  class_means <- matrix(1:4, nrow = n_class, ncol = n_gene)

  # This works ONLY because all classes are of the same size
  # Otherwise the mean of the means != overall mean
  overall_means <- colSums(class_means) / n_class
  psd <- matrix(1, nrow = 1, ncol = n_gene)
  stats <- calculate_class_stats(classes, class_means, overall_means, psd)
  expect_equal(
    stats[1, 1],
    (class_means[1, 1] - overall_means[1]) /
      (sqrt(1 / n_rep - 1 / length(classes)) * psd[1])
  )
})
