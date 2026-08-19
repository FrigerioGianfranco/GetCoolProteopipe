#' Change the name of samples.
#'
#' Change the name of some samples.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#' @param old_names character containing the name of samples you want to modify.
#' @param new_names character containing the new name of samples.
#' @param old_new_names_table alternatively, you can bass here a table with old names in the first column and new names in the second column. If you pass an argument here, this will be considered instead of old_names and new_names.
#'
#' @return the GCPlist with the sample names updated in the sampleINFO, and intensities data frames.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist01 <- GCP_ChangeSampleNames(GCPlist = GCPlist00,
#'                                    old_names = GCPlist00$sampleINFO$Sample,
#'                                    new_names = c("S_1", "S_2", "S_3", "S_4", "S_5", "V_1", "V_2", "V_3", "V_4", "V_5"))
#'
#' }
#'
#'
#' @export
GCP_ChangeSampleNames <- function(GCPlist, old_names = pull(GCPlist$sampleINFO, 1), new_names = paste0(old_names, "_updated"), old_new_names_table = NULL) {

  checkGCPlist(GCPlist)


  if (is.null(old_new_names_table)) {

    if (any(is.na(old_names))) {stop("old_names must not contain missing values")}

    if (length(old_names) != length(new_names)) {stop("old_names and new_names must have the same length")}


    old_names <- as.character(old_names)
    new_names <- as.character(new_names)

    if (!all(old_names %in% pull(GCPlist$sampleINFO, 1))) {
      stop(paste0("Names passed to old_names must be names of samples. The following sample names you passed are not present:\n",
                  paste0(old_names[which(!old_names %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
    }

    if (any(is.na(new_names))) {
      new_names[which(is.na(new_names))] <- old_names[which(is.na(new_names))]
    }


    name_changing_table <- tibble(old = old_names, new = new_names)

  } else {

    old_new_names_table <- as.data.frame(old_new_names_table)

    if (ncol(old_new_names_table) != 2) {stop("old_new_names_table must contain exactly two columns")}

    if (any(is.na(pull(old_new_names_table, 1)))) {stop("the first column of old_new_names_table must not contain missing values")}

    name_changing_table <- tibble(old = as.character(pull(old_new_names_table, 1)),
                                  new = as.character(pull(old_new_names_table, 2)))

    if (!all(name_changing_table$old %in% pull(GCPlist$sampleINFO, 1))) {
      stop(paste0("Names passed to the first column of old_new_names_table must be names of samples. The following sample names you passed are not present:\n",
                  paste0(name_changing_table$old[which(!name_changing_table$old %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
    }

    if (any(is.na(name_changing_table$new))) {
      name_changing_table$new[which(is.na(name_changing_table$new))] <- name_changing_table$old[which(is.na(name_changing_table$new))]
    }
  }

  GCPoutput <- GCPlist

  for (i in 1:nrow(name_changing_table)) {

    GCPoutput$sampleINFO[which(pull(GCPoutput$sampleINFO, 1) == name_changing_table$old[i]), 1] <- name_changing_table$new[i]
    colnames(GCPoutput$intensities)[which(colnames(GCPoutput$intensities) == name_changing_table$old[i])] <- name_changing_table$new[i]

  }


  cat("\n")
  cat("______\n")
  cat("The name of samples have been updated in this way:\n\n")

  data_frame_to_print <- data.frame(ORIGINAL_NAMES = pull(GCPlist$sampleINFO, 1),
                                    NAMES_UPDATED = pull(GCPoutput$sampleINFO, 1),
                                    COMMENT = as.character(NA))
  data_frame_to_print[,"COMMENT"] <- ifelse(data_frame_to_print$ORIGINAL_NAMES != data_frame_to_print$NAMES_UPDATED, "  Changed!", "  Not changed")
  data_frame_to_print[,"ORIGINAL_NAMES"] <- ifelse(data_frame_to_print$COMMENT == "  Not changed", "", data_frame_to_print$ORIGINAL_NAMES)

  print(data_frame_to_print)

  cat("______\n")

  return(GCPoutput)

}


