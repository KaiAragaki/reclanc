## code to prepare `synthetic_expression` dataset goes here
library(DESeq2)

dds <- makeExampleDESeqDataSet(
  n = 100,
  m = 12,
  betaSD = 1
)

expression <- dds |>
  normTransform() |>
  assay()

synthetic_expression <- list(expression = expression, classes = dds$condition)

usethis::use_data(synthetic_expression, overwrite = TRUE)
