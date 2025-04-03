#' Replacing zero with NA
#'
#' In the quant_raw and the quantLFQ data intensity, it replaces all zeros with missing values.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#'
#' @return a GCPlist list with zero replaced with NAs in the data intensity tables.
#'
#' @export
GCP_ReplaceZerowithNA <- function(GCPlist) {

  checkGCPlist(GCPlist)


  GCPoutput <- GCPlist

  cat("\n")
  cat("______\n")
  cat("The number of zeros replaced with NAs is\n")

  counting_tot_raw <- 0
  total_cells_raw <- nrow(GCPoutput$quant_raw)*length(colnames(GCPoutput$quant_raw)[-which(colnames(GCPoutput$quant_raw) == "protid")])

  for (a in colnames(GCPoutput$quant_raw)[-which(colnames(GCPoutput$quant_raw) == "protid")]) {

    counting_sa <- 0

    if (length(which(pull(GCPoutput$quant_raw, a) == 0))>0) {
      for (i in which(pull(GCPoutput$quant_raw, a) == 0)) {
        GCPoutput$quant_raw[i, a] <- NA
        counting_sa <- counting_sa+1
        counting_tot_raw <- counting_tot_raw+1
      }
    }
  }


  cat(paste0(" - ", counting_tot_raw, " out of ", total_cells_raw, " (", round(counting_tot_raw/total_cells_raw*100, digits = 1), "%) in the raw table;\n"))


  counting_tot_LFQ <- 0
  total_cells_LFQ <- nrow(GCPoutput$quant_LFQ)*length(colnames(GCPoutput$quant_LFQ)[-which(colnames(GCPoutput$quant_LFQ) == "protid")])

  for (a in colnames(GCPoutput$quant_LFQ)[-which(colnames(GCPoutput$quant_LFQ) == "protid")]) {

    counting_sa <- 0

    if (length(which(pull(GCPoutput$quant_LFQ, a) == 0))>0) {
      for (i in which(pull(GCPoutput$quant_LFQ, a) == 0)) {
        GCPoutput$quant_LFQ[i, a] <- NA
        counting_sa <- counting_sa+1
        counting_tot_LFQ <- counting_tot_LFQ+1
      }
    }
  }

  cat(paste0(" - ", counting_tot_LFQ, " out of ", total_cells_LFQ, " (", round(counting_tot_LFQ/total_cells_LFQ*100, digits = 1), "%) in the LFQ table.\n"))

  cat("______\n")

  return(GCPoutput)

}
