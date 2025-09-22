#' ratio PTM on protein intensities
#'
#' It takes the each valid feature intensity of the quant_raw table generated from PTM and it divide it by the relative raw intensity of the protein; the match is based on the Accession.
#'
#' @param GCPlistPTM a list created with the ImportPTMs function.
#' @param GCPlistProteins a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". If "raw", the denominator of the ratio will be the quant_raw table of GCPlistProteins; if "LFQ", it will be the quantLFQ table.
#'
#' @return a GCPlist list with the calculated ratios of intensities in the quant_raw table.
#'
#' @export
GCP_ratioPTMprotein <- function(GCPlistPTM, GCPlistProteins, raw_or_LFQ = c("lfq", "raw")) {

  checkGCPlist(GCPlistPTM)
  checkGCPlist(GCPlistProteins)


  if (!"Accession"%in%colnames(GCPlistPTM$proteinINFO)) {stop("The column Accession is not present in the proteinINFO table of GCPlistPTM")}
  if (!"Accession"%in%colnames(GCPlistProteins$proteinINFO)) {stop("The column Accession is not present in the proteinINFO table of GCPlistProteins")}

  if (any(duplicated(GCPlistProteins$proteinINFO$Accession))) {
    stop(paste0("The coulumn Accession of the proteinINFO table of GCPlistProteins contains duplicates! In particular, the followings are duplicated:\n",
                paste0(GCPlistProteins$proteinINFO$Accession[which(duplicated(GCPlistProteins$proteinINFO$Accession))], collapse = "\n")))
  }

  if (!identical(tolower(raw_or_LFQ), c("lfq", "raw"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("lfq", "raw"))

  if (raw_or_LFQ == "lfq") {
    cat("\n -- LFQ data are used for proteins as denominator --\n\n")
  } else if (raw_or_LFQ == "raw") {
    cat("\n -- raw data are used for proteins as denominator --\n\n")
  }


  all_protidPTM <- GCPlistPTM$quant_raw$protid
  valid_protidPTM <- character()
  negative_protidPTM <- character()

  for (i in 1:length(all_protidPTM)) {
    this_prot_values <- as.numeric(GCPlistPTM$quant_raw[which(GCPlistPTM$quant_raw$protid==all_protidPTM[i]), which(colnames(GCPlistPTM$quant_raw)!="protid")])
    this_prot_values_noNA <- this_prot_values[which(!is.na(this_prot_values))]
    this_prot_values_noNA_noZERO <- this_prot_values_noNA[which(this_prot_values_noNA!=0)]

    if (any(this_prot_values_noNA_noZERO<0)) {
      negative_protidPTM <- c(negative_protidPTM, all_protidPTM[i])
    }


    if (length(this_prot_values_noNA_noZERO)>0) {
      valid_protidPTM <- c(valid_protidPTM, all_protidPTM[i])
    }
  }

  if (length(negative_protidPTM)>0) {
    stop(paste0("The following protid in PTM contain negative values!\n",
                paste0(negative_protidPTM, collapse = "\n")))
  }

  if (length(valid_protidPTM)<1) {
    stop("There are no valid intensities in the PTM to use!!")
  }

  GCPoutput_nr <- GCPlistPTM

  tibble_validPTM <- tibble(protid = valid_protidPTM,
                            Accession = GCPlistPTM$proteinINFO$Accession[which(GCPlistPTM$proteinINFO$protid%in%valid_protidPTM)])



  if (!all(tibble_validPTM$Accession %in% GCPlistProteins$proteinINFO$Accession)) {

    tibble_validPTM_not_found <- filter(tibble_validPTM,
                                        !Accession %in% GCPlistProteins$proteinINFO$Accession)

    tibble_validPTM_found <- filter(tibble_validPTM,
                                    Accession %in% GCPlistProteins$proteinINFO$Accession)

    GCPoutput_nr$quant_raw <- filter(GCPoutput_nr$quant_raw,
                        protid %in% tibble_validPTM_found$protid)
    GCPoutput_nr$quant_LFQ <- filter(GCPoutput_nr$quant_LFQ,
                                  protid %in% tibble_validPTM_found$protid)

  } else {

    tibble_validPTM_not_found <- tibble_validPTM[0,]

    tibble_validPTM_found <- tibble_validPTM
  }

  tibble_validPTM_found_paired <- add_column(tibble_validPTM_found, protidProteins = as.character(NA))

  for (i in 1:nrow(tibble_validPTM_found_paired)) {
    tibble_validPTM_found_paired[i, "protidProteins"] <- GCPlistProteins$proteinINFO$protid[which(GCPlistProteins$proteinINFO$Accession == tibble_validPTM_found_paired$Accession[i])]
  }

  if (raw_or_LFQ == "lfq") {
    protein_table_to_use <- GCPlistProteins$quant_LFQ[which(GCPlistProteins$quant_LFQ$protid %in% tibble_validPTM_found_paired$protidProteins),]
  } else if (raw_or_LFQ == "raw") {
    protein_table_to_use <- GCPlistProteins$quant_raw[which(GCPlistProteins$quant_raw$protid %in% tibble_validPTM_found_paired$protidProteins),]
  } else {
    stop('raw_or_LFQ must be "raw" or "LFQ"')
  }



  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(is.na(.))))) {
    stop("There are missing values in the protein intensities! Before using this function, at least use GCP_NAimputation!")
  }

  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(.<0)))) {
    stop("There are negative values among the protein intensities!")
  }

  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(.==0)))) {
    stop("There are values equal to zero in the protein intensities! Before using this function, at least use GCP_ReplaceZerowithNA and then GCP_NAimputation!")
  }


  if (!all(colnames(GCPoutput_nr$quant_raw)[which(colnames(GCPoutput_nr$quant_raw)!="protid")] %in% colnames(protein_table_to_use)[which(colnames(protein_table_to_use)!="protid")])) {
    stop(paste0("The following samples are not present in the protein table:\n",
                paste0(colnames(GCPoutput_nr$quant_raw)[which(colnames(GCPoutput_nr$quant_raw)!="protid")][which(!colnames(GCPoutput_nr$quant_raw)[which(colnames(GCPoutput_nr$quant_raw)!="protid")] %in% colnames(protein_table_to_use)[which(colnames(protein_table_to_use)!="protid")])],
                       collapse = "\n")))
  }


  if (!all(tibble_validPTM$Accession %in% GCPlistProteins$proteinINFO$Accession)) {
    cat(paste0("\n\n  ", nrow(tibble_validPTM_not_found), " out of ", nrow(tibble_validPTM), " has been removed as there is no correspondence with any Accession in the protein table:\n"))
    print(as.data.frame(tibble_validPTM_not_found))
  } else {
    cat("\n\n  Found a correspondence for all the Accension of the protein table!")
  }
  cat("\n")



  GCPoutput <- GCPoutput_nr

  for (o in 1:nrow(GCPoutput_nr$quant_raw)) {

    this_protidPTM <- GCPoutput_nr$quant_raw$protid[o]
    this_protidProteins <- tibble_validPTM_found_paired$protidProteins[which(tibble_validPTM_found_paired$protid==this_protidPTM)]

    for (a in colnames(GCPoutput_nr$quant_raw)[which(colnames(GCPoutput_nr$quant_raw)!="protid")]) {

      if (!is.na(pull(GCPoutput_nr$quant_raw, a)[o])) {
        GCPoutput$quant_raw[o, a] <- pull(GCPoutput_nr$quant_raw, a)[o]/pull(protein_table_to_use, a)[which(protein_table_to_use$protid==this_protidProteins)]
      }

    }
  }

  return(GCPoutput)
}

