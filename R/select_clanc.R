#' @param dk a matrix of statistics, where each column is a class and each row
#'   is a gene
#' @param active number of active genes to consider. Can either be a single
#'   number (in which case the same number of active genes will be used for each
#'   class) or a vector of numbers of length equal to the number of classes,
#'   with the order corresponding to the levels of the classes factor.
clanc_select <- function(dk, active) {
  active <- validate_active(dk, active)
  active <- data.frame(
    class = levels(dk$class),
    active = active
  )
  all <- dplyr::left_join(dk, active, by = "class")
  all$reserved <- FALSE
  all$n_win <- 0
  selection_recurse(all) |>
    dplyr::filter(gets) |>
    dplyr::select(class, gene_id)
}

selection_recurse <- function(df) {
  if (all(df$reserved) || all(df$n_win >= df$active)) return(NULL)
  df <- df |>
    filter(!reserved, n_win < active) |>
    arrange(gene_id, rank, desc(abs_stat)) |>
    mutate(n = 1:n(), win = n == 1, .by = "gene_id") |>
    arrange(class, rank) |>
    mutate(
      n_win = n_win + cumsum(win),
      gets = win & (n_win <= active),
      .by = "class"
    ) |>
    mutate(reserved = any(gets), .by = "gene_id") |>
    mutate(n_win = max(n_win), .by = "class")
  rbind(df, selection_round(df))
}

validate_active <- function(dk, active) {
  # NOTE will this cause a problem if a give fold doesn't have a level of a
  # class?
  n_classes <- length(levels(dk$class))

  if (!length(active) %in% c(1, n_classes)) stop("Invalid `length(active)`")

  df <- aggregate(no_tie ~ classes, dk, sum)
  df <- df[order(df$classes), ] # Order rows by class factor...just in case.
  df$active <- active
  df$n_genes <- length(unique(dk$gene_id))

  if (any(df$active > df$n_genes)) {
    warning("# active genes > # genes available in at least one class")
    message("Using all genes for those classes.")
  }
  if (any(df$active > df$no_tie)) {
    warning("# active genes > # uniquely ranked genes in at least one class")
    message("Reducing active gene number for this class to max possible.")
  }
  mapply(min, df$no_ties, df$n_genes, df$active)
}
