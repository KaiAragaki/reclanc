length_1_and_char <- function(x) {
  length(x) == 1 && is.character(x)
}

# TRUE if the user seems like they're asking for a column in colData
# AND it's there

# FALSE if it doesn't seem like the user is asking for something in colData

# ERROR if user seems like they're asking for column in colData and it isn't
# there
spec_in_cd <- function(x, cd_names) {
  if (!length_1_and_char(x)) {
    return(FALSE)
  }
  if (!x %in% cd_names) {
    cli::cli_abort("{x} is not a column name in {.code colData(x)}")
  }
  TRUE
}
