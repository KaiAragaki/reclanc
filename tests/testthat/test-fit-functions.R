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
