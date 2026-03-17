#' Set the name_column_groups option.
#'
#' It sets the option GetCoolProteopipe.name_column_groups, so it will be consistenly applied to all the functions that have the name_column_groups argument.
#'
#' @param name_column_groups NULL or a character of length 1. The option GetCoolProteopipe.name_column_groups will be set to it.
#'
#' @return Nothing. It only sets the GetCoolProteopipe.name_column_groups option.
#'
#'
#' @examples
#' \dontrun{
#'
#' set_name_column_groups("Condition")
#'
#' }
#'
#'
#' @export
set_name_column_groups <- function(name_column_groups = NULL) {

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups) != 1) {stop("name_column_groups must be NULL or a charcter of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a charcter of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a charcter of length 1, not a missing value")}

    cat(paste0("\n --- the name_column_groups is now set to be '", name_column_groups, "' ---\n\n"))

  } else {
    cat(paste0("\n --- the name_column_groups is now set to be NULL ---\n\n"))
  }

  options(GetCoolProteopipe.name_column_groups = name_column_groups)
  invisible(getOption("GetCoolProteopipe.name_column_groups"))
}
