mark_active_genes <- function(expression, active) {
  active <- validate_active(expression, active)
  active <- data.frame(
    class = levels(expression$class),
    active = active
  )
  dplyr::left_join(expression, active, by = "class") |>
    dplyr::mutate(reserved = FALSE, n_win = 0) |>
    selection_recurse() |>
    dplyr::filter(.data$gets)
}

# The meat of the select function
# A recursive function.
#
# In the base case, all genes are 'reserved' (they have been claimed by a class)
# or all classes have claimed as many genes as they want (as set by 'active').
#
# For the base case, return NULL
#
# Otherwise...
#
# Get rid of any reserved genes, and get rid of any classes who are already full
# up on active genes
#
# All classes hold a tournament for each gene. The class that has the highest
# 'personal rank' for that gene gets 1st place. If there are ties, the absolute
# value of the statistic breaks it. The class that gets first place is the
# winner, all others are losers.
#
# After this is done for all genes, the genes are arranged by class and rank.
#
# Classes only take as many genes as they need, and other their top picks. If a
# class wins more genes than it needs, it will recuse itself from further rounds
# and release the unneeded genes.
#
# Any genes that were needed and taken are marked as reserved.
#
# This result is recursively rbind'd to calls of itself.
selection_recurse <- function(expression) {
  if (all(expression$reserved) || all(expression$n_win >= expression$active)) {
    return(NULL)
  }

  expression <- expression |>
    dplyr::filter(
      !.data$reserved,
      .data$n_win < .data$active
    ) |>
    dplyr::arrange(.data$gene_id, .data$rank, dplyr::desc(.data$abs_stat)) |>
    dplyr::mutate(
      n = seq_len(dplyr::n()),
      win = .data$n == 1,
      .by = "gene_id"
    ) |>
    dplyr::arrange(.data$class, .data$rank) |>
    dplyr::mutate(
      n_win = .data$n_win + cumsum(.data$win),
      gets = .data$win & (.data$n_win <= .data$active),
      .by = "class"
    ) |>
    dplyr::mutate(reserved = any(.data$gets), .by = "gene_id") |>
    dplyr::mutate(n_win = max(.data$n_win), .by = "class")
  rbind(expression, selection_recurse(expression))
}

validate_active <- function(expression, active) {
  # NOTE will this cause a problem if a give fold doesn't have a level of a
  # class?
  n_classes <- length(levels(expression$class))

  if (!length(active) %in% c(1, n_classes)) stop("Invalid `length(active)`")

  df <- aggregate(no_tie ~ class, expression, sum)
  df <- df[order(df$class), ] # Order rows by class factor...just in case.
  df$active <- active
  df$n_genes <- length(unique(expression$gene_id))

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
