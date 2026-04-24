#' Set the log_base option.
#'
#' It sets the option GetCoolProteopipe.log_base, so it will be consistently applied to all the functions that have the base or log_base argument.
#'
#' @param log_base NULL or a number. The option GetCoolProteopipe.log_base will be set to it.
#'
#' @return Nothing. It only sets the GetCoolProteopipe.log_base option.
#'
#'
#' @examples
#' \dontrun{
#'
#' set_log_base(2)
#'
#' }
#'
#'
#' @export
set_log_base <- function(log_base = 2) {

  if (!is.null(log_base)) {
    if (length(log_base) != 1) {stop("log_base must be NULL or a numeric of length 1")}
    if (!is.numeric(log_base)) {stop("log_base must be NULL or a numeric of length 1")}
    if (is.na(log_base)) {stop("log_base must be NULL or a numeric of length 1, not a missing value")}

    cat(paste0("\n --- the log_base option is now set to be ", log_base, " ---\n\n"))

  } else {
    cat(paste0("\n --- the log_base option is now set to be NULL ---\n\n"))
  }

  options(GetCoolProteopipe.log_base = log_base)
  invisible(getOption("GetCoolProteopipe.log_base"))
}
