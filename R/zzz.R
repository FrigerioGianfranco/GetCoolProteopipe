# file zzz.R

#' @import tidyverse
#' @import GetFeatistics

.onLoad <- function(libname, pkgname) {

  op <- options()
  op_pkg <- list(
    GetCoolProteopipe.name_column_groups = NULL,
    GetCoolProteopipe.col_pal = NULL,
    GetCoolProteopipe.log_base = 2
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

  if (!is.null(getOption("GetCoolProteopipe.name_column_groups"))) {
    cat(paste0("\n --- the name_column_groups is set to be '", getOption("GetCoolProteopipe.name_column_groups"), "' ---\n\n"))
  }

  if (!is.null(getOption("GetCoolProteopipe.col_pal"))) {
    if (is.null(names(getOption("GetCoolProteopipe.col_pal")))) {
      cat(paste0('\n --- the col_pal option is set to be  c("', paste0(getOption("GetCoolProteopipe.col_pal"), collapse = '", "'), '") ---\n\n'))
    } else {
      cat('\n --- the col_pal option is set to be  c(')
      for (i in seq(length(getOption("GetCoolProteopipe.col_pal")))) {
        if (names(getOption("GetCoolProteopipe.col_pal"))[i]!="") {
          cat(paste0(names(getOption("GetCoolProteopipe.col_pal"))[i], ' = "'))
        } else {
          cat('"')
        }
        cat(paste0(getOption("GetCoolProteopipe.col_pal")[i], '"'))
        if (i != length(getOption("GetCoolProteopipe.col_pal"))) {
          cat(', ')
        }
      }
      cat(') ---\n\n')
    }
  }

  if (!is.null(getOption("GetCoolProteopipe.log_base"))) {
    if (getOption("GetCoolProteopipe.log_base")!=2) {
      cat(paste0("\n --- the log_base option is set to be ", getOption("GetCoolProteopipe.log_base"), " ---\n\n"))
    }
  }
}
