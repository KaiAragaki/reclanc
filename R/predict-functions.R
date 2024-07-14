filter_predictors <- function(fit, predictors) {
  # make sure rownames exist
  # Check number of genes in centroids that are in dataset. if no rows after
  # filter, provide info of what the centroid names look like (example).
  merge(fit$centroids, predictors, by.x = "gene", by.y = "row.names")
}

calc_dist <- function(all) {
  user_col_ids <- seq(from = 7, to = ncol(all))
  user_data <- all[, user_col_ids, drop = FALSE]

  vv <- all$pooled_sd^2
  first_term <- (all$expression / all$pooled_sd)^2 |>
    stats::aggregate(by = list(all$class), FUN = sum)
  second_term <- (-2 * (user_data * all$expression / vv)) |>
    stats::aggregate(by = list(all$class), FUN = sum)

  priors <- unique(data.frame(class = all$class, prior = all$prior))
  priors$log_priors <- 2 * log(priors$prior)

  dists <- sweep(second_term[-1], 1, first_term[[2]], FUN = "+") |>
    sweep(1, priors$log_priors)

  classes <- factor(levels(priors$class), levels = levels(priors$class))

  data.frame(class = classes, dists)
}

calc_cors <- function(all, method) {
  splits <- split(all, ~class)
  do.call(
    hardhat::spruce_numeric_multiple,
    lapply(splits, calc_cor, method = method)
  )
}

calc_cor <- function(split, method) {
  cor <- stats::cor(split$expression, split[7:ncol(split)], method = method)
  rownames(cor) <- unique(split$class)
  as.numeric(cor)
}

predict_class <- function(dists) {
  dists[["class"]][apply(dists[2:ncol(dists)], 2, which.min)]
}

predict_dist_multi <- function(dists) {
  rownames(dists) <- dists$class
  dists <- dists[-1]
  do.call(hardhat::spruce_numeric_multiple, as.list(data.frame(t(dists))))
}

# dists: col = sample, row = class, value = dist

wrangle_data <- function(data) {
  if (inherits(data, "SummarizedExperiment")) return(wrangle_data_se(data))
  if (inherits(data, "data.frame")) return(wrangle_data_df(data))
  if (inherits(data, "matrix")) return(wrangle_data_matrix(data))
  if (inherits(data, "ExpressionSet")) return(wrangle_data_es(data))
  cli::cli_abort("Unable to predict object of type {class(data)}")
}

wrangle_data_se <- function(data) {
  # Should assay number be an optional arg?
  data <- SummarizedExperiment::assay(data)
  t(data)
}

wrangle_data_df <- function(data) {
  data
}

wrangle_data_matrix <- function(data) {
  data
}

wrangle_data_es <- function(data) {
  t(Biobase::exprs(data))
}


# Often, not all of the predictors will be present
# This is typical and worth working around rather than erroring
custom_forge <- function(new_data, blueprint) {
  predictors <- colnames(blueprint$ptypes$predictors)
  contains <- intersect(predictors, colnames(new_data))

  if (length(contains) == 0) {
    cli::cli_abort("No genes in centroids detected in `new_data`")
  }

  u_pred <- length(unique(predictors))
  u_nd <- length(contains)
  pct_found <- round(100 * u_nd/ u_pred)

  cli::cli_inform(
    "{u_nd}/{u_pred} ({pct_found}%) genes in centroids found in `new_data`"
  )

  new_data <- new_data[, which(colnames(new_data) %in% contains)]
  new_data <- tibble::as_tibble(new_data)
  list(predictors = new_data, outcomes = NULL, extras = NULL)
}
