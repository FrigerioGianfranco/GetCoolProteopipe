#' Replacing zero with NA
#'
#' In the intensities table, it replaces all zeros with missing values.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#'
#' @return a GCPlist list with zero replaced with NAs in the data intensity table.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist04 <- GCP_ReplaceZerowithNA(GCPlist03)
#'
#' }
#'
#'
#' @export
GCP_ReplaceZerowithNA <- function(GCPlist) {

  checkGCPlist(GCPlist)

  GCPoutput <- GCPlist

  cat("\n")
  cat("______\n")
  cat("The number of zeros replaced with NAs is\n")

  counting_tot <- 0
  total_cells <- nrow(GCPoutput$intensities)*length(colnames(GCPoutput$intensities)[-which(colnames(GCPoutput$intensities) == "protid")])

  for (a in colnames(GCPoutput$intensities)[-which(colnames(GCPoutput$intensities) == "protid")]) {

    counting_sa <- 0

    if (length(which(pull(GCPoutput$intensities, a) == 0))>0) {
      for (i in which(pull(GCPoutput$intensities, a) == 0)) {
        GCPoutput$intensities[i, a] <- NA
        counting_sa <- counting_sa+1
        counting_tot <- counting_tot+1
      }
    }
  }

  cat(paste0(" - ", counting_tot, " out of ", total_cells, " (", round(counting_tot/total_cells*100, digits = 1), "%).\n"))

  cat("______\n")

  return(GCPoutput)
}
