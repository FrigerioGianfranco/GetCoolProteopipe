#' Set the use of raw or LFQ.
#'
#' It sets the option GetCoolProteopipe.raw_or_LFQ, so it will be consistenly applied to all the functions that have the raw_or_LFQ argument.
#'
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The option GetCoolProteopipe.raw_or_LFQ will be set to it.
#'
#' @return Nothing. It only sets the GetCoolProteopipe.raw_or_LFQ option.

#' @export
set_raw_or_LFQ <- function(raw_or_LFQ = c("lfq", "raw")) {

  if (!identical(tolower(raw_or_LFQ), c("lfq", "raw"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("lfq", "raw"))

  if (raw_or_LFQ == "lfq") {
    cat("\n --- the option is now set to use LFQ data ---\n\n")
  } else if (raw_or_LFQ == "raw") {
    cat("\n --- the option is now set to use raw data ---\n\n")
  }

  options(GetCoolProteopipe.raw_or_LFQ = raw_or_LFQ)
  invisible(getOption("GetCoolProteopipe.raw_or_LFQ"))
}
