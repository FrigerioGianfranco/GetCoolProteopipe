# file zzz.R

#' @import tidyverse
#' @import GetFeatistics

.onLoad <- function(libname, pkgname) {
  packageStartupMessage("\nGetCoolProteopipe v", packageVersion("GetCoolProteopipe"), "\n")

  library(tidyverse)
  library(GetFeatistics)
}
