new_clanc <- function(centroids, blueprint) {
  hardhat::new_model(
    centroids = centroids,
    blueprint = blueprint,
    class = "clanc"
  )
}

validate_centroids <- function(centroids, class_data) {
  stopifnot(is.matrix(centroids))
  if (any(duplicated(rownames(centroids)))) {
    stop("centroids must not have duplicate rownames")
  }
}

validate_class_data <- function(centroids, class_data) {
  stopifnot(is.data.frame(class_data))
  stopifnot(ncol(centroids) == nrow(class_data))
  stopifnot(all.equal(colnames(class_data), c("class", "active", "priors")))
  stopifnot(is.factor(class_data$class))
  stopifnot(is.numeric(class_data$active))
  stopifnot(is.numeric(class_data$priors))
}

validate_standard_deviations <- function(centroids, standard_deviations) {
  stopifnot(is.vector(standard_deviations))
}
