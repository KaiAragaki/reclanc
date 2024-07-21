clanc_fit <- function(expression, class_data, classes) {
  class_means <- collapse::fmean(expression, classes, na.rm = TRUE)
  overall_means <- colMeans(expression, na.rm = TRUE)

  pooled_sds <- calculate_pooled_sd(
    expression,
    class_data,
    classes,
    class_means
  )

  class_stats <- calculate_class_stats(
    classes, class_means, overall_means, pooled_sds
  )
  abs_stats <- abs(class_stats)

  ranks <- t(apply(abs_stats, 1, \(x) rank(-x, ties.method = "min")))
  ties <- t(apply(ranks, 1, duplicated))
  selected <- select_genes(abs_stats, ranks, ties, class_data)
  centroids <- create_centroids(selected, class_means, overall_means)

  pooled_sds <- pooled_sds[match(rownames(centroids), names(pooled_sds))]

  tall_centroids <- matrix(unlist(centroids), ncol = 1) |> as.data.frame()
  tall_centroids$gene <- rownames(centroids)
  tall_centroids$class <- rep(colnames(centroids), each = nrow(centroids))
  tall_centroids <- tall_centroids[, c(3, 2, 1)]
  w_sd <- merge(tall_centroids, pooled_sds, by.x = "gene", by.y = "row.names")
  out <- merge(w_sd, class_data, by = "class")
  out$class <- factor(out$class, levels = levels(class_data$class))
  colnames(out) <- c(
    "class", "gene", "expression", "pooled_sd", "active", "prior"
  )

  list(centroids = out)
}

#' @importFrom collapse %c*%
calculate_pooled_sd <- function(expression, class_data, classes, class_means) {
  class_means <- as.matrix(class_means)
  df <- nrow(expression) - nrow(class_data)

  mat <- matrix(0, nrow = 1, ncol = ncol(expression))
  for (i in seq_along(levels(classes))) {
    class_exp <- expression[which(classes == levels(classes)[i]), ]
    fswept <- collapse::TRA(class_exp, class_means[i, , drop = FALSE])
    fsquared <- fswept %c*% fswept
    fsummed <- collapse::fsum(fsquared)
    mat <- mat + fsummed
  }

  out <- sqrt(colSums(mat / df))
  names(out) <- colnames(expression)
  out
}

#' @importFrom collapse %c/% %r-% %r/%
calculate_class_stats <- function(classes,
                                  class_means,
                                  overall_means,
                                  class_pooled_sds) {
  mks <- sqrt(1 / table(classes) - 1 / length(classes))
  class_means <- as.matrix(class_means)
  fsweep1 <- class_means %r-% overall_means
  fsweep2 <- fsweep1 %r/% class_pooled_sds
  fsweep3 <- fsweep2 %c/% as.data.frame(mks)$Freq
  fsweep3
}

select_genes <- function(abs_stats, ranks, ties, class_data) {
  df <- data.frame(
    gene = colnames(abs_stats),
    abs = matrix(t(abs_stats), ncol = 1),
    rank = matrix(t(ranks), ncol = 1),
    tie = matrix(t(ties), ncol = 1),
    class = rep(class_data$class, each = ncol(abs_stats)),
    reserved = FALSE,
    n_win = 0
  )

  # class_data is where `active` lives
  df <- merge(df, class_data, by = "class") |>
    selection_recurse()

  df <- df[df$gets, ]

  data.frame(
    class = factor(df$class, levels = levels(class_data$class)),
    gene = df$gene
  )
}



#' @importFrom rlang .data
selection_recurse <- function(df) {
  if (all(df$reserved) || all(df$n_win >= df$active)) {
    return(NULL)
  }

  df <- df |>
    dplyr::filter(!.data$reserved, .data$n_win < .data$active) |>
    dplyr::arrange(.data$gene, .data$rank, dplyr::desc(.data$abs)) |>
    dplyr::mutate(
      n = seq_len(dplyr::n()),
      win = .data$n == 1,
      .by = "gene"
    ) |>
    dplyr::arrange(.data$class, .data$rank) |>
    dplyr::mutate(
      n_win = .data$n_win + cumsum(.data$win),
      gets = .data$win & (.data$n_win <= .data$active),
      .by = "class"
    ) |>
    dplyr::mutate(reserved = any(.data$gets), .by = "gene") |>
    dplyr::mutate(n_win = max(.data$n_win), .by = "class")
  rbind(df, selection_recurse(df))
}

create_centroids <- function(selected, class_means, overall_means) {
  lvls <- levels(selected$class)
  n_levels <- length(lvls)
  class <- class_means[, match(selected$gene, colnames(class_means))]
  overall <- overall_means[match(selected$gene, names(overall_means))]
  mm <- matrix(rep(overall, n_levels))

  update_col <- function(class_means, winner, overall_means) {
    idx <- which(winner == lvls)
    new <- rep(overall_means, length(class_means))
    new[idx] <- class_means[idx]
    new
  }

  out <- mapply(update_col, class, selected$class, overall)
  out <- t(out)
  colnames(out) <- lvls
  out
}

# Predictors to use are only 'discovered' *after* fitting
# Therefore the default way of dealing with predictors isn't sufficient
make_new_ptypes <- function(fit, processed) {
  genes <- unique(fit$centroids$gene)
  # make a dummy tibble w proper dims
  preds <- matrix(1.0, nrow = 1, ncol = length(genes))
  colnames(preds) <- genes
  preds <- tibble::as_tibble(preds)
  preds <- preds[-1, ]
  list(predictors = preds, outcomes = processed$blueprint$ptypes$outcomes)
}
