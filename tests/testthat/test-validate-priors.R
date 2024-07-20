test_that("priors processed correctly", {
  processed <- list(outcomes = data.frame(class = factor(c("a", "b", "a"))))
  expect_error(process_priors(processed, data.frame()))
  expect_equal(process_priors(processed, "equal"), c(0.5, 0.5))
  expect_equal(process_priors(processed, "class"), c(2 / 3, 1 / 3))
  expect_equal(process_priors(processed, c(0.2, 0.8)), c(0.2, 0.8))
})

test_that("char priors processed correctly", {
  expect_equal(
    process_priors_char("equal", factor(c("a", "b", "c", "a"))),
    c(1 / 3, 1 / 3, 1 / 3)
  )
  expect_equal(
    process_priors_char("class", factor(c("a", "b", "c", "a"))),
    c(1 / 2, 1 / 4, 1 / 4)
  )
  expect_equal(
    process_priors_char(
      "class",
      factor(c("a", "b", "c", "a"), levels = c("c", "b", "a"))
    ),
    c(1 / 4, 1 / 4, 1 / 2)
  )
  expect_equal(process_priors_char("equal", factor("a")), 1)
  expect_equal(process_priors_char("equal", factor(c("a", "a"))), 1)
  expect_error(process_priors_char(c("equal", "equal"), factor("a")))
})

test_that("num priors processed correctly", {
  expect_equal(process_priors_num(1, factor("a")), 1)
  expect_warning(process_priors_num(0.9, factor("a")))
  expect_error(process_priors_num(-1, factor("a")))
  expect_error(process_priors_num(2, factor("a")))
  expect_error(process_priors_num(c(-1, 2), factor(c("a", "b"))))
  expect_warning(process_priors_num(c(0.5, 0.4), factor(c("a", "b"))))
  expect_error(process_priors_num(c(0.5, 0.6), factor(c("a", "b"))))
  expect_error(process_priors_num(1, factor(c("a", "b"))))
  expect_no_error(process_priors_num(1.000001, factor("a")))
  expect_error(process_priors_num(1, factor(c("a", "b"))))
})
