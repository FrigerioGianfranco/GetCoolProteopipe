#' ratio PTM on protein intensities
#'
#' It takes each valid feature intensity from PTM and it divides it by the relative intensity of the protein; the match is based on the Accession.
#'
#' @param GCPlistPTM a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#' @param GCPlistProteins a list created with the ImportOutputMaxQuant function.
#' @param are_log_transf logical. If the intensities are log-transformed, specify here as TRUE, so the subtraction will be performed instead of the ratio.
#'
#' @return a GCPlist list with the calculated ratios of intensities in the intensities table.
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' GCPlistPTM09 <- GCP_ratioPTMprotein(GCPlistPTM = GCPlistPTM08,
#'                                     GCPlistProteins = GCPlist14,
#'                                     are_log_transf = TRUE)
#'
#'
#' }
#'
#'
#' @export
GCP_ratioPTMprotein <- function(GCPlistPTM, GCPlistProteins, are_log_transf = TRUE) {

  checkGCPlist(GCPlistPTM)
  checkGCPlist(GCPlistProteins)

  if (length(are_log_transf)!=1) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (!is.logical(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (is.na(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}


  if (!"Accession"%in%colnames(GCPlistPTM$proteinINFO)) {stop("The column Accession is not present in the proteinINFO table of GCPlistPTM")}
  if (!"Accession"%in%colnames(GCPlistProteins$proteinINFO)) {stop("The column Accession is not present in the proteinINFO table of GCPlistProteins")}

  if (any(duplicated(GCPlistProteins$proteinINFO$Accession))) {
    stop(paste0("The coulumn Accession of the proteinINFO table of GCPlistProteins contains duplicates! In particular, the followings are duplicated:\n",
                paste0(GCPlistProteins$proteinINFO$Accession[which(duplicated(GCPlistProteins$proteinINFO$Accession))], collapse = "\n")))
  }

  all_protidPTM <- GCPlistPTM$intensities$protid
  valid_protidPTM <- character()
  negative_protidPTM <- character()

  for (i in 1:length(all_protidPTM)) {
    this_prot_values <- as.numeric(GCPlistPTM$intensities[which(GCPlistPTM$intensities$protid==all_protidPTM[i]), which(colnames(GCPlistPTM$intensities)!="protid")])
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



  Protein_proteinINFO_filtered <- GCPlistProteins$proteinINFO[which(GCPlistProteins$proteinINFO$protid%in%GCPlistProteins$intensities$protid),]



  GCPoutput_nr <- GCPlistPTM

  tibble_validPTM <- tibble(protid = valid_protidPTM,
                            Accession = GCPlistPTM$proteinINFO$Accession[which(GCPlistPTM$proteinINFO$protid%in%valid_protidPTM)])



  if (!all(tibble_validPTM$Accession %in% Protein_proteinINFO_filtered$Accession)) {

    tibble_validPTM_not_found <- filter(tibble_validPTM,
                                        !Accession %in% Protein_proteinINFO_filtered$Accession)

    tibble_validPTM_found <- filter(tibble_validPTM,
                                    Accession %in% Protein_proteinINFO_filtered$Accession)

    GCPoutput_nr$intensities <- filter(GCPoutput_nr$intensities,
                                       protid %in% tibble_validPTM_found$protid)

  } else {

    tibble_validPTM_not_found <- tibble_validPTM[0,]

    tibble_validPTM_found <- tibble_validPTM
  }

  tibble_validPTM_found_paired <- add_column(tibble_validPTM_found, protidProteins = as.character(NA))

  for (i in 1:nrow(tibble_validPTM_found_paired)) {
    tibble_validPTM_found_paired[i, "protidProteins"] <- Protein_proteinINFO_filtered$protid[which(Protein_proteinINFO_filtered$Accession == tibble_validPTM_found_paired$Accession[i])]
  }


  protein_table_to_use <- GCPlistProteins$intensities[which(GCPlistProteins$intensities$protid %in% tibble_validPTM_found_paired$protidProteins),]




  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(is.na(.))))) {
    stop("There are missing values in the protein intensities! Before using this function, at least use GCP_NAimputation!")
  }

  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(.<0)))) {
    stop("There are negative values among the protein intensities!")
  }

  if (any(map_lgl(protein_table_to_use[,which(colnames(protein_table_to_use)!="protid")], ~ any(.==0)))) {
    stop("There are values equal to zero in the protein intensities! Before using this function, at least use GCP_ReplaceZerowithNA and then GCP_NAimputation!")
  }


  if (!all(colnames(GCPoutput_nr$intensities)[which(colnames(GCPoutput_nr$intensities)!="protid")] %in% colnames(protein_table_to_use)[which(colnames(protein_table_to_use)!="protid")])) {
    stop(paste0("The following samples are not present in the protein table:\n",
                paste0(colnames(GCPoutput_nr$intensities)[which(colnames(GCPoutput_nr$intensities)!="protid")][which(!colnames(GCPoutput_nr$intensities)[which(colnames(GCPoutput_nr$intensities)!="protid")] %in% colnames(protein_table_to_use)[which(colnames(protein_table_to_use)!="protid")])],
                       collapse = "\n")))
  }


  if (!all(tibble_validPTM$Accession %in% Protein_proteinINFO_filtered$Accession)) {
    cat(paste0("\n\n  ", nrow(tibble_validPTM_not_found), " out of ", nrow(tibble_validPTM), " has been removed as there is no correspondence with any Accession in the protein table:\n"))
    print(as.data.frame(tibble_validPTM_not_found))
  } else {
    cat("\n\n  Found a correspondence for all the Accension of the protein table!")
  }
  cat("\n")



  GCPoutput <- GCPoutput_nr

  for (o in 1:nrow(GCPoutput_nr$intensities)) {

    this_protidPTM <- GCPoutput_nr$intensities$protid[o]
    this_protidProteins <- tibble_validPTM_found_paired$protidProteins[which(tibble_validPTM_found_paired$protid==this_protidPTM)]

    for (a in colnames(GCPoutput_nr$intensities)[which(colnames(GCPoutput_nr$intensities)!="protid")]) {

      if (!is.na(pull(GCPoutput_nr$intensities, a)[o])) {
        if (!are_log_transf) {
          GCPoutput$intensities[o, a] <- pull(GCPoutput_nr$intensities, a)[o]/pull(protein_table_to_use, a)[which(protein_table_to_use$protid==this_protidProteins)]
        } else {
          GCPoutput$intensities[o, a] <- pull(GCPoutput_nr$intensities, a)[o]-pull(protein_table_to_use, a)[which(protein_table_to_use$protid==this_protidProteins)]
        }
      }
    }
  }

  return(GCPoutput)
}
