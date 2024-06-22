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
    aggregate(by = list(all$class), FUN = sum)
  second_term <- (-2 * (user_data * all$expression / vv)) |>
    aggregate(by = list(all$class), FUN = sum)

  priors <- unique(data.frame(class = all$class, prior = all$prior))
  priors$log_priors <- 2 * log(priors$prior)

  dists <- sweep(second_term[-1], 1, first_term[[2]], FUN = "+") |>
    sweep(1, priors$log_priors)

  classes <- factor(levels(priors$class), levels = levels(priors$class))

  data.frame(class = classes, dists)
}

predict_class <- function(dists) {
  dists[["class"]][apply(dists[2:ncol(dists)], 2, which.min)]
}
