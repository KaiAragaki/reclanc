test_that("catch unsupported types", {
  expect_error(clanc(numeric()))
})

test_that("SummarizedExperiment catch bad args", {
  library(SummarizedExperiment)
  assay <- matrix(1:4, nrow = 2)
  cd <- DataFrame(
    a = c("one", "two"),
    active = c(3, 4),
    prior = c(0.5, 0.5)
  )
  rd <- DataFrame(gene = c("g1", "g2"))
  se <- SummarizedExperiment(assay, rowData = rd, colData = cd)
  expect_error(clanc(se, classes = "b"))
  expect_error(clanc(se, active = "active", classes = c("one", "two")))
  expect_error(clanc(se, priors = "prior", classes = c("one", "two")))

  # Too many attributes per class:
  assay <- matrix(1:6, nrow = 2)
  cd <- DataFrame(
    a = c("one", "two", "one"),
    active = c(3, 4, 2),
    prior = c(0.5, 0.5, 0.5)
  )
  rd <- DataFrame(gene = c("g1", "g2"))
  se <- SummarizedExperiment(assay, rowData = rd, colData = cd)
  expect_error(
    clanc(se, classes = "a", active = "active", priors = "prior"),
    "More than one prior"
  )
})

test_that("ExpressionSet catch bad args", {
  library(Biobase)
  assay <- matrix(1:4, nrow = 2)
  pdata <- AnnotatedDataFrame(
    data.frame(a = c("one", "two"), active = c(3, 4), prior = c(0.5, 0.5))
  )
  fdata <- AnnotatedDataFrame(data.frame(gene = c("g1", "g2")))
  es <- ExpressionSet(assay, pdata, fdata)
  expect_error(clanc(es, classes = "b"))
  expect_error(clanc(es, active = "active", classes = c("one", "two")))
  expect_error(clanc(es, priors = "prior", classes = c("one", "two")))

  # Too many attributes per class:
  assay <- matrix(1:6, nrow = 2)
  pdata <- AnnotatedDataFrame(
    data.frame(
      a = c("one", "two", "one"),
      active = c(3, 4, 2),
      prior = c(0.5, 0.5, 0.5)
    )
  )
  fdata <- AnnotatedDataFrame(data.frame(gene = c("g1", "g2")))
  es <- ExpressionSet(assay, pdata, fdata)
  expect_error(
    clanc(es, classes = "a", active = "active", priors = "prior"),
    "More than one prior"
  )
})

test_that("fit works", {
  data <- synthetic_expression

  fit_mat <- clanc(data$expression, data$classes, active = 2)

  fit_df <- clanc(as.data.frame(data$expression), data$classes, active = 2)

  se <- SummarizedExperiment::SummarizedExperiment(
    data$expression,
    colData = S4Vectors::DataFrame(class = data$classes)
  )

  fit_se <- clanc(se, classes = "class", active = 2, assay = 1)

  pdata <- Biobase::AnnotatedDataFrame(
    data.frame(class = data$class,row.names = colnames(data$expression))
  )

  es <- Biobase::ExpressionSet(data$expression, pdata)

  fit_es <- clanc(es, classes = "class", active = 2)

  form <- cbind(class = data$classes, as.data.frame(t(data$expression)))

  form_fit <- clanc(class ~ ., form, active = 2)

  recipe_fit <- parsnip::discrim_linear() |>
    parsnip::set_engine("clanc", active = 2) |>
    parsnip::fit(class ~ ., data = form)

  expect_true(
    all(fit_mat$centroids$gene %in% c("gene13", "gene41", "gene52", "gene74"))
  )

  expect_true(
    all(
      sapply(
        list(fit_df, fit_se, fit_es, form_fit, recipe_fit$fit),
        identical,
        fit_mat
      )
    )
  )
})
