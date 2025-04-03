#' Remove samples with all zero intensities
#'
#' Starting from a GCPlist, it removes proteins which intensity is equal to zero in all the samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param modality one of the following: "raw_only", "LFQ_only", "either", "both". If "raw_only", it filters rows that are all zero in the raw table; if "LFQ_only", it filters rows that are all zero in the LFQ table; if "either", it filters rows that are all zero in either the raw or the LFQ table; if "both", it filters rows that are all zero in both the raw and LFQ tables.
#'
#' @return the GCPlist in which each table has potentially a reduced number of rows.
#'
#' @export
GCP_RemoveAllZero <- function(GCPlist, modality = c("raw_only", "LFQ_only", "either", "both")) {

  checkGCPlist(GCPlist)

  if (!identical(tolower(modality), c("raw_only", "lfq_only", "either", "both"))) {
    if (length(modality) != 1) {stop('modality must be one of "raw_only", "LFQ_only", "either", "both"')}
    if (is.na(modality)) {stop('modality must be one of "raw_only", "LFQ_only", "either", "both"')}
  }
  modality <- tolower(modality)
  modality <- match.arg(modality, c("raw_only", "lfq_only", "either", "both"))


  protid_total <- GCPlist$proteinINFO$protid
  protid_zero_in_raw <- character()
  protid_zero_in_LFQ <- character()


  for (i in 1:length(protid_total)) {
    if (all(as.numeric(GCPlist$quant_raw[i, -1])==0)) {
      protid_zero_in_raw <- c(protid_zero_in_raw, GCPlist$quant_raw$protid[i])
    }

    if (all(as.numeric(GCPlist$quant_LFQ[i, -1])==0)) {
      protid_zero_in_LFQ <- c(protid_zero_in_LFQ, GCPlist$quant_LFQ$protid[i])
    }
  }

  protid_zero_in_either <- unique(c(protid_zero_in_raw, protid_zero_in_LFQ))
  protid_zero_in_both <- protid_zero_in_either[which(protid_zero_in_either%in%protid_zero_in_raw & protid_zero_in_either%in%protid_zero_in_LFQ)]
  protid_zero_uniquely_in_raw <- protid_zero_in_either[which(protid_zero_in_either%in%protid_zero_in_raw & !protid_zero_in_either%in%protid_zero_in_LFQ)]
  protid_zero_uniquely_in_LFQ <- protid_zero_in_either[which(!protid_zero_in_either%in%protid_zero_in_raw & protid_zero_in_either%in%protid_zero_in_LFQ)]

  cat("\n")
  cat("___\n")
  cat(paste0("There are ", length(protid_zero_in_raw), " rows with all zero in the raw intensity table.\n"))
  cat(paste0("There are ", length(protid_zero_in_LFQ), " rows with all zero in the LFQ intensity table.\n"))
  cat(paste0("(", length(protid_zero_in_both), " from both tables, ", length(protid_zero_uniquely_in_raw), " uniquely from raw, and ", length(protid_zero_uniquely_in_LFQ), " uniquely from LFQ)\n"))
  cat("\n")

  if (tolower(modality) == "either") {

    GCPoutput <- GCPlist
    GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_zero_in_either),]
    GCPoutput$quant_raw <- GCPoutput$quant_raw[which(!GCPoutput$quant_raw$protid%in%protid_zero_in_either),]
    GCPoutput$quant_LFQ <- GCPoutput$quant_LFQ[which(!GCPoutput$quant_LFQ$protid%in%protid_zero_in_either),]

    cat(paste0('...since you chose the modality "either", ', length(protid_zero_in_either), ' proteins were removed.\n'))
    cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept\n.'))


  } else if (tolower(modality) == "raw_only") {

    GCPoutput <- GCPlist
    GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_zero_in_raw),]
    GCPoutput$quant_raw <- GCPoutput$quant_raw[which(!GCPoutput$quant_raw$protid%in%protid_zero_in_raw),]
    GCPoutput$quant_LFQ <- GCPoutput$quant_LFQ[which(!GCPoutput$quant_LFQ$protid%in%protid_zero_in_raw),]

    cat(paste0('...since you chose the modality "raw_only", ', length(protid_zero_in_raw), ' proteins were removed.\n'))
    cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept\n.'))

  } else if (tolower(modality) == "lfq_only") {

    GCPoutput <- GCPlist
    GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_zero_in_LFQ),]
    GCPoutput$quant_raw <- GCPoutput$quant_raw[which(!GCPoutput$quant_raw$protid%in%protid_zero_in_LFQ),]
    GCPoutput$quant_LFQ <- GCPoutput$quant_LFQ[which(!GCPoutput$quant_LFQ$protid%in%protid_zero_in_LFQ),]

    cat(paste0('...since you chose the modality "LFQ_only", ', length(protid_zero_in_LFQ), ' proteins were removed.\n'))
    cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept\n.'))

  } else if (tolower(modality) == "both") {

    GCPoutput <- GCPlist
    GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_zero_in_both),]
    GCPoutput$quant_raw <- GCPoutput$quant_raw[which(!GCPoutput$quant_raw$protid%in%protid_zero_in_both),]
    GCPoutput$quant_LFQ <- GCPoutput$quant_LFQ[which(!GCPoutput$quant_LFQ$protid%in%protid_zero_in_both),]

    cat(paste0('...since you chose the modality "both", ', length(protid_zero_in_both), ' proteins were removed.\n'))
    cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept.\n'))

  }

  cat("___\n")

  return(GCPoutput)

}
