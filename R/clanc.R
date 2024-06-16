#' @param expression a matrix of expression data, where columns are samples and
#'   rows are genes.
#' @param classes a factor vector of classes for each column in `expression`
#' @param priors either a numeric vector of length 1, ncol(expression), "equal",
#'   or "class". If "equal", equal probabilities will be used for each class. If
#'   "class", the proportions of each class in the training data will be used as
#'   the prior.
#' @param active number of genes to use. If a single number, the same number of
#'   genes will be used for every class. Alternatively, a different number can
#'   be given for each class by supplying a vector.

clanc <- function(expression,
                  classes,
                  priors = "equal",
                  active = 10) {
  stopifnot(
    is.factor(classes),
    is.matrix(expression),
    ncol(expression) == length(classes)
  )


}

clanc.ExpressionSet <- function() {

}

clanc.SummarizedExperiment <- function() {

}

clanc.default <- function() {

}

# Metadata approach - class, prior, active table?
# The full classes vector would need to be supplied as well.
