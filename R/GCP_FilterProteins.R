#' Filtering proteins
#'
#' It filters the quant_raw and quant_LFQ data frames considering certain criteria.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param operation character. Operation(s) to be applied considering a column of proteinINFO (for example: "ttest_PvaluesFDR < 0.05"). All the TRUE from the operation will be kept.
#'
#' @return The GCPlist with a potentially reduced number of rows in the quant_raw and quant_LFQ tables.
#'
#' @export
GCP_FilterProteins <- function(GCPlist, operation) {

  checkGCPlist(GCPlist)

  if (!is.character(operation)) {stop("operation must be a character vector")}
  if (length(operation) < 1) {stop("operation must contain at least something to do!")}
  if (any(is.na(operation))) {stop("operation must not contain NAs")}


  SUMMARY_DATAFRAME <- data.frame(quant_raw = nrow(GCPlist$quant_raw), quant_LFQ = nrow(GCPlist$quant_LFQ))
  row.names(SUMMARY_DATAFRAME) <- "Initially, the number of proteins were "


  GCPout <- GCPlist

  for (op in operation) {

    proteinINFO_fil <- filter(GCPout$proteinINFO, eval(parse(text = op)))

    GCPout$quant_raw <- filter(GCPout$quant_raw, protid %in% proteinINFO_fil$protid)
    GCPout$quant_LFQ <- filter(GCPout$quant_LFQ, protid %in% proteinINFO_fil$protid)

    SUMMARY_DATAFRAME <- rbind(SUMMARY_DATAFRAME,
                               data.frame(quant_raw = nrow(GCPout$quant_raw), quant_LFQ = nrow(GCPout$quant_LFQ)))
    row.names(SUMMARY_DATAFRAME)[nrow(SUMMARY_DATAFRAME)] <- paste0("after ", op, " ")
  }


  cat("\n")

  print(SUMMARY_DATAFRAME)

  cat("\n")


  return(GCPout)
}
