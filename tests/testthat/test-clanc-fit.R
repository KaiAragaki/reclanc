test_that("catch unsupported types", {
  expect_error(clanc(numeric()))
})

test_that("catch bad args", {
  library(SummarizedExperiment)
  assay <- matrix(1:4, nrow = 2)
  cd <- DataFrame(a = c("one", "two"), active = c(3, 4), prior = c(0.5, 0.5))
  rd <- DataFrame(gene = c("g1", "g2"))
  se <- SummarizedExperiment(assay, rowData = rd, colData = cd)
  expect_error(clanc(se, classes = "b"))
  expect_error(clanc(se, active = "active", classes = c("one", "two")))
  expect_error(clanc(se, priors = "prior", classes = c("one", "two")))
  # Too many attributes per class:
  assay <- matrix(1:6, nrow = 2)
  cd <- DataFrame(a = c("one", "two", "one"), active = c(3, 4, 2), prior = c(0.5, 0.5, 0.5))
  rd <- DataFrame(gene = c("g1", "g2"))
  se <- SummarizedExperiment(assay, rowData = rd, colData = cd)
  expect_error(
    clanc(se, classes = "a", active = "active", priors = "prior"),
    "More than one prior"
  )
})