# file zzz.R

#' @import tidyverse
#' @import GetFeatistics

.onLoad <- function(libname, pkgname) {

  op <- options()
  op_pkg <- list(
    GetCoolProteopipe.raw_or_LFQ = "lfq"
  )
  toset <- !(names(op_pkg) %in% names(op))
  if (any(toset)) {options(op_pkg[toset])}
  invisible()

  library(tidyverse)
  suppressPackageStartupMessages(
    library(GetFeatistics, quietly = TRUE, warn.conflicts = FALSE)
  )
}


.onAttach <- function(libname, pkgname) {
  packageStartupMessage("\nGetCoolProteopipe v", packageVersion("GetCoolProteopipe"), "\n")

  if (getOption("GetCoolProteopipe.raw_or_LFQ") == "lfq") {
    cat("\n --- the option is set to use LFQ data ---\n\n")
  } else if (getOption("GetCoolProteopipe.raw_or_LFQ") == "raw") {
    cat("\n --- the option is set to use raw data ---\n\n")
  }
}
