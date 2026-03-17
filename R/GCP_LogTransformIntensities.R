#' Log-transform intensities
#'
#' It log-transforms all the intensities of the intensities table.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param base numeric. The base of the logarithm to use.
#'
#' @return a GCPlist list with log-transformed data intensity table.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist06 <- GCP_LogTransformIntensities(GCPlist = GCPlist05,
#'                                          base = 2)
#'
#' }
#'
#'
#' @export
GCP_LogTransformIntensities <- function(GCPlist, base = exp(1)) {

  checkGCPlist(GCPlist)

  if (length(base)!=1) {stop("base must be a numeric of length 1")}
  if (!is.numeric(base)) {stop("base must be a numeric of length 1")}
  if (is.na(base)) {stop("base must be a numeric of length 1, not a missing value")}

  if (any(map_lgl(GCPlist$intensities[, -which(colnames(GCPlist$intensities) == "protid")], \(x) any(x[which(!is.na(x))] == 0)))) {
    cat("\n")
    cat("\n")
    cat("There are some zeros in the the intensities! Are you sure you don't want to use the ReplaceZerowithNA function first?!")
    cat("\n")
  }

  if (any(map_lgl(GCPlist$intensities[, -which(colnames(GCPlist$intensities) == "protid")], \(x) any(x[which(!is.na(x))] < 0)))) {
    cat("\n")
    cat("\n")
    cat("There are negative numbers in the the intensities! There is no logarithm of a negative number, so keep in mind that NaNs will be introduced...!!")
    cat("\n")
  }


  GCPoutput <- GCPlist

  for (a in colnames(GCPoutput$intensities)[-which(colnames(GCPoutput$intensities) == "protid")]) {

    for (i in 1:length(pull(GCPoutput$intensities, a))) {
      GCPoutput$intensities[i, a] <- log(pull(GCPoutput$intensities, a)[i], base = base)
    }
  }

  if (base == exp(1)) {
    cat("\n______\nAll intensities have been transformed using the natural logarithm.\n______\n")
  } else {
    cat(paste0("\n______\nAll intensities have been transformed using the logarithm base ", base, ".\n______\n"))
  }

  return(GCPoutput)
}

