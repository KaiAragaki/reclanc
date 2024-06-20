filter_predictors <- function(fit, predictors) {
  # make sure rownames exist
  # Check number of genes in centroids that are in dataset. if no rows after
  # filter, provide info of what the centroid names look like (example).
  merge(fit$centroids, predictors, by.x = "gene", by.y = "row.names")
}

calc_dist <- function(all) {
  user_col_ids <- seq(from = 7, to = ncol(all))
  user_data <- all[, user_col_ids, drop = FALSE]
  zs <- user_data |>
    sweep(1, all$expression) |>
    sweep(1, all$pooled_sd, FUN = "/")

  priors <- unique(data.frame(class = all$class, prior = all$prior))
  priors$log_priors <- 2 * log(priors$prior)
  dist <- aggregate(zs^2, by = list(all$class), FUN = sum)
  merged <- merge(priors, dist, by.x = "class", by.y = "Group.1")
  user_ids <- seq(from = 4, to = ncol(merged))
  merged[, user_ids] <- sweep(
    merged[, user_ids, drop = FALSE],
    1,
    merged$log_priors
  )
  merged
}

predict_class <- function(dists) {
  dists[["class"]][apply(dists[4:ncol(dists)], 2, which.min)]
}
