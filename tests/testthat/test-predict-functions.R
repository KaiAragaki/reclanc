test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("predict works", {
  data <- synthetic_expression
  fit_mat <- clanc(data$expression, data$classes, active = 2)
  preds <- predict(fit_mat, data$expression, "class", format = "tall")$.pred_class
  preds_cor <- predict(fit_mat, data$expression, "numeric", format = "tall", method = "pearson")
  preds_cor <- c("A", "B")[apply(preds_cor, 1, which.max)] |> factor()
  truth <- factor(rep(c("A", "B"), each = 6))
  expect_equal(preds, truth)
  expect_equal(preds_cor, truth)
})
