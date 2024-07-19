#' Calculate centroids from expression data with ClaNC
#'
#' @param x Depending on the context:
#'
#'   * A __data frame__ of expression.
#'   * A __matrix__ of expression.
#'   * A __recipe__ specifying a set of preprocessing steps
#'     created from [recipes::recipe()].
#'   * An __ExpressionSet__.
#'   * A __SummarizedExperiment__ with `assay` containing expression.
#'
#' Expression should be library-size corrected, but not scaled.
#'
#' If supplying a __data frame__, __matrix__, __ExpressionSet__,
#' __SummarizedExperiment__, the rows should represent genes, and the columns
#' should represent samples (as is standard for expression data). The column
#' names should be sample IDs, while the row names should be gene IDs.
#'
#' If a __recipe__ is provided, the data should have genes as
#' columns (to match the formula provided to the recipe.)
#'
#' @param classes When `x` is a __data frame__ or __matrix__, `class` contains
#'   class labels with the form of either:
#'
#'   * A __data frame__ with 1 factor column
#'   * A factor __vector__.
#'
#' When `x` is an __ExpressionSet__ or __SummarizedExperiment__, `class` is the
#' name of the column in `pData(x)` or `colData(x)` that contains classes as a
#' factor.
#'
#' @param data When a __recipe__ or __formula__ is used, `data` is specified as:
#'
#'   * A __data frame__ containing both expression and classes, where columns
#'   are the genes or class, and rows are the samples.
#'
#' @param assay When a __SummarizedExperiment__ is used, the index or name of
#'   the assay
#'
#' @param formula A formula specifying the classes on the left-hand side,
#' and the predictor terms on the right-hand side.
#'
#' @param active Either a single number or a numeric vector equal to the length
#'   of the number of unique class labels. Represents the number class-specific
#'   genes that should be selected for a centroid. Note that different numbers
#'   of genes can be selected for each class. See details.
#'
#' When `x` is an __ExpressionSet__ or __SummarizedExperiment__, `active` can
#' additionally by the name of the column in `pData(x)` or `colData(x)` that
#' contains the numeric vector
#'
#' @param priors Can take a variety of values:
#'
#'   * "equal" - each class has an equal prior
#'   * "class" - each class has a prior equal to its frequency in the training
#'   set
#'   * A numeric vector with length equal to number of classes
#'
#' When `x` is an __ExpressionSet__ or __SummarizedExperiment__, `active` can
#' additionally by the name of the column in `pData(x)` or `colData(x)` that
#' contains the numeric vector
#'
#' @param ... Not currently used, but required for extensibility.
#'
#' @details The original description of ClaNC can be found
#'   [here](10.1093/bioinformatics/bti681)
#'
#' While `active` sets the number of class-specific genes, each centroid will
#' have more than that number of genes. To explain by way of example, if `active
#' = 5` and there are 3 classes, each centroid will have 15 genes, with 5 of
#' those genes being particular to a given class. If these genes are 'active' in
#' that class, their values will be the mean of the class. If the genes are not
#' active in that given class, their values will be the overall expression of
#' the given gene across all classes.
#'
#' @return
#'
#' A `clanc` object.
#'
#' @examples
#'
#' expression_matrix <- synthetic_expression$expression
#' head(expression_matrix)
#' classes <- synthetic_expression$classes
#' classes
#'
#' # data.frame/tibble/matrix interface:
#'
#' clanc(expression_matrix, classes = classes, active = 5, priors = "equal")
#'
#
#' # Formula interface:
#'
#' # Data must have class included as a column
#' # Genes must be *columns* and samples must be *rows*
#' # Hence the data transposition.
#' for_formula <- data.frame(class = classes, t(expression_matrix))
#'
#' clanc(class ~ ., for_formula, active = 5, priors = "equal")
#'
#'
#' # Recipes interface:
#'
#' rec <- recipes::recipe(class ~ ., data = for_formula)
#'
#' clanc(rec, for_formula, active = 5, priors = "equal")
#'
#' # SummarizedExperiment interface:
#' se <- SummarizedExperiment::SummarizedExperiment(
#'   expression_matrix,
#'   colData = data.frame(
#'     class = classes,
#'     active = 5,
#'     prior = c(0.5, 0.5)
#'   )
#' )
#'
#' clanc(se, classes = "class", active = "active", priors = "equal")
#'
#' # ExpressionSet interface:
#' adf <- data.frame(
#'   row.names = colnames(expression_matrix),
#'   class = classes
#' ) |>
#'   Biobase::AnnotatedDataFrame()
#'
#' es <- Biobase::ExpressionSet(expression_matrix, adf)
#
#' clanc(es, classes = "class", active = 5, priors = 0.5)
#'
#' @export
clanc <- function(x, ...) {
  UseMethod("clanc")
}

#' @export
#' @rdname clanc
clanc.default <- function(x, ...) {
  cli::cli_abort("`clanc()` is not defined for a {class(x)[1]}.")
}

# XY method - data frame

#' @export
#' @rdname clanc
clanc.data.frame <- function(x,
                             classes,
                             active,
                             priors = "equal",
                             ...) {
  x <- t(x)
  processed <- hardhat::mold(x, classes)
  clanc_bridge(processed, active, priors, ...)
}

# XY method - matrix

#' @export
#' @rdname clanc
clanc.matrix <- function(x,
                         classes,
                         active,
                         priors = "equal",
                         ...) {
  x <- t(x)
  processed <- hardhat::mold(x, classes)
  clanc_bridge(processed, active, priors, ...)
}

# XY method - SummarizedExperiment

#' @export
#' @rdname clanc
clanc.SummarizedExperiment <- function(x,
                                       classes,
                                       active,
                                       priors = "equal",
                                       assay = 1,
                                       ...) {
  rlang::check_installed(
    "SummarizedExperiment", reason = "to use `clanc.SummarizedExperiment()`"
  )
  expression <- SummarizedExperiment::assay(x, assay)
  expression <- t(expression)
  cd_names <- colnames(SummarizedExperiment::colData(x))

  if (!spec_in_cd(classes, cd_names) &&
        (spec_in_cd(priors, cd_names) || spec_in_cd(active, cd_names))) {
    cli::cli_abort(
      "`classes` must be specified as a column name in colData if `active` or `priors` are." #nolint
    )
  }

  if (spec_in_cd(classes, cd_names))
    metadata <- data.frame(class = x[[classes]])

  if ((!priors %in% c("equal", "class")) && spec_in_cd(priors, cd_names))
    metadata <- cbind(metadata, prior = x[[priors]])

  if (spec_in_cd(active, cd_names))
    metadata <- cbind(metadata, active = x[[active]])

  if (spec_in_cd(classes, cd_names) &&
        nrow(unique(metadata)) > length(unique(metadata[[classes]]))) {
    cli::cli_abort(c(
      "More than one prior and/or active for each class",
      "i" = "Each class must have exactly 1 prior and active",
      "i" = "To test across many `active`, use `tune::tune_grid`"
    ))
  } else {
    classes <- metadata$class
    metadata <- unique(metadata)
    if ((!priors %in% c("equal", "class")) && spec_in_cd(priors, cd_names)) {
      priors <- metadata$priors
    }
    if (spec_in_cd(active, cd_names)) active <- metadata$active
  }

  processed <- hardhat::mold(expression, classes)
  clanc_bridge(processed, active, priors, ...)
}

# XY method - ExpressionSet

#' @export
#' @rdname clanc
clanc.ExpressionSet <- function(x,
                                classes,
                                active,
                                priors = "equal",
                                ...) {
  rlang::check_installed(
    "Biobase", reason = "to use `clanc.ExpressionSet()`"
  )
  expression <- t(Biobase::exprs(x))
  pd_names <- colnames(Biobase::pData(x))

  if (!spec_in_cd(classes, pd_names) &&
        (spec_in_cd(priors, pd_names) || spec_in_cd(active, pd_names))) {
    cli::cli_abort(
      "`classes` must be specified as a column name in pData if `active` or `priors` are." #nolint
    )
  }
  if (spec_in_cd(classes, pd_names))
    metadata <- data.frame(class = x[[classes]])

  if ((!priors %in% c("equal", "class")) && spec_in_cd(priors, pd_names))
    metadata <- cbind(metadata, prior = x[[priors]])

  if (spec_in_cd(active, pd_names))
    metadata <- cbind(metadata, active = x[[active]])

  if (spec_in_cd(classes, pd_names) &&
        nrow(unique(metadata)) > length(unique(metadata$class))) {
    cli::cli_abort(
      "More than one prior and/or active for each class",
      "i" = "Each class must have exactly 1 prior and active",
      "i" = "To test across many `active`, use `tune::tune_grid`"
    )
  } else {
    metadata <- unique(metadata)
    classes <- x[[classes]]
    if ((!priors %in% c("equal", "class")) && spec_in_cd(priors, pd_names)) {
      priors <- metadata$priors
    }
    if (spec_in_cd(active, pd_names)) active <- metadata$active
  }

  processed <- hardhat::mold(expression, classes)
  clanc_bridge(processed, active, priors, ...)
}

# Formula method

#' @export
#' @rdname clanc
clanc.formula <- function(formula,
                          data,
                          active,
                          priors = "equal",
                          ...) {
  # Oftentimes the formula will be incredibly wide.
  # R hates this.
  #
  # We get around this using some recent functions from recipes and pretend a
  # data.frame was supplied
  args <- form2args(formula, data)
  data <- args$data
  pred <- data[, which(colnames(data) %in% args$predictors)]
  outcomes <- data[, which(colnames(data) %in% args$outcomes)]

  clanc.data.frame(t(pred), outcomes, active, priors, ...)
}

# Recipe method

#' @export
#' @rdname clanc
clanc.recipe <- function(x,
                         data,
                         active,
                         priors = "equal",
                         ...) {
  processed <- hardhat::mold(x, data)
  clanc_bridge(processed, active, priors, ...)
}

# ------------------------------------------------------------------------------
# Bridge

clanc_bridge <- function(processed, active, priors, ...) {

  # NOTE What happens if classes supplied upstream are not the proper length?
  # That is, the length of the vector is shorter than the length of the data?
  # (not necessarily n_levels deficiency)
  active <- process_active(processed, active, priors, ...)
  priors <- process_priors(processed, active, priors, ...)
  validate_classes(processed, active, priors)

  expression <- processed$predictors
  # hardhat 'mold' returns diff colnames for recipe and all others
  classes <- processed$outcomes[[1]]
  class_data <- data.frame(
    class = sort(unique(classes)),
    active = active,
    priors = priors
  )
  fit <- clanc_impl(expression, class_data, classes)

  # Create a new blueprint using only active genes
  bp <- hardhat::new_xy_blueprint(
    intercept = FALSE,
    allow_novel_levels = FALSE,
    composition = "tibble",
    ptypes = make_new_ptypes(fit, processed)
  )

  new_clanc(
    centroids = fit$centroids,
    blueprint = bp
  )
}

# ------------------------------------------------------------------------------
# Implementation

clanc_impl <- function(expression, class_data, classes) {
  fit <- clanc_fit(expression, class_data, classes)
  list(centroids = fit$centroids)
}
