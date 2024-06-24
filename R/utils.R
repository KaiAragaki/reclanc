get_var_info_from_form <- function(se, formula, var_n) {
  formula_names <- all.vars(formula)
  name <- formula_names[var_n]
  data <- SummarizedExperiment::colData(se)[[name]]
  list(name = name, data = data)
}
