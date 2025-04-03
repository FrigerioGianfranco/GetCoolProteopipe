#' Log-transform intensities
#'
#' In the quant_raw and the quantLFQ data intensity, it log-transforms all the intensities.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param base numeric. The base of the logarithm to use.
#'
#' @return a GCPlist list with log-transformed data intensity tables.
#'
#' @export
GCP_LogTransformIntensities <- function(GCPlist, base = exp(1)) {

  checkGCPlist(GCPlist)

  if (length(base)!=1) {stop("base must be a numeric of length 1")}
  if (!is.numeric(base)) {stop("base must be a numeric of length 1")}
  if (is.na(base)) {stop("base must be a numeric of length 1, not a missing value")}

  if (any(map_lgl(GCPlist$quant_raw[, -which(colnames(GCPlist$quant_raw) == "protid")], \(x) any(x[which(!is.na(x))] == 0))) | any(map_lgl(GCPlist$quant_LFQ[, -which(colnames(GCPlist$quant_LFQ) == "protid")], \(x) any(x[which(!is.na(x))] == 0)))) {
    cat("\n")
    cat("\n")
    cat("There are some zeros in the the intensities! Are you sure you don't want to use the ReplaceZerowithNA function first?!")
    cat("\n")
  }

  if (any(map_lgl(GCPlist$quant_raw[, -which(colnames(GCPlist$quant_raw) == "protid")], \(x) any(x[which(!is.na(x))] < 0))) | any(map_lgl(GCPlist$quant_LFQ[, -which(colnames(GCPlist$quant_LFQ) == "protid")], \(x) any(x[which(!is.na(x))] < 0)))) {
    cat("\n")
    cat("\n")
    cat("There are negative numbers in the the intensities! There is no logarithm of a negative number, so keep in mind that NaNs will be introduced...!!")
    cat("\n")
  }



  GCPoutput <- GCPlist

  for (a in colnames(GCPoutput$quant_raw)[-which(colnames(GCPoutput$quant_raw) == "protid")]) {

    for (i in 1:length(pull(GCPoutput$quant_raw, a))) {
      GCPoutput$quant_raw[i, a] <- log(pull(GCPoutput$quant_raw, a)[i], base = base)

    }
  }


  for (a in colnames(GCPoutput$quant_LFQ)[-which(colnames(GCPoutput$quant_LFQ) == "protid")]) {

    for (i in 1:length(pull(GCPoutput$quant_LFQ, a))) {
      GCPoutput$quant_LFQ[i, a] <- log(pull(GCPoutput$quant_LFQ, a)[i], base = base)

    }
  }

  return(GCPoutput)

}

