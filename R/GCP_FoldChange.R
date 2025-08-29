#' Performing a Fold Change analysis on data
#'
#' It performs a Fold Change analyses on the proteins intensities. Please, be aware that the Fold Change analysis should be performed only on positive data! Indeed, all protein intensities should be positive, or unreliable results will be generated!
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function. IMPORTANT: All protein intensities should be positive, or unreliable results will be generated!
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The Fold Change analysis will be performed only in the specified data intensities.
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups. Since this is a Fold Change analysis, there must be exactly two groups.
#' @param control_group NULL or a character of length 1. The name of the group to be considered as control group, which will be considered as denominator in the ratio of the Fold Change analysis. If NULL, the first level of the factor will be considered as the control group.
#' @param paired logical. If FALSE it performs FC on mean of the two groups. If TRUE it performs FC for each pair and then compute the mean.
#' @param are_log_transf logical. If the protein intensities are already log-transformed, specify here as TRUE, so the subtraction will be performed instead of the ratio.
#' @param log_base numeric of length 1. Specify here the base of the logarithm to calculate the logFC or, if are_log_transf is TRUE; the base of the logarithm that were used to transform the data.
#'
#' @return The GCPlist with the results of the Fold Change analysis added to the proteinINFO data frame.
#'
#' @export
GCP_FoldChange <- function(GCPlist, raw_or_LFQ = c("lfq", "raw"), name_column_groups, control_group = NULL,
                           paired = FALSE, are_log_transf = TRUE, log_base = 2) {

  checkGCPlist(GCPlist)

  if (!identical(tolower(raw_or_LFQ), c("lfq", "raw"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("lfq", "raw"))

  if (raw_or_LFQ == "lfq") {
    cat("\n -- LFQ data are used --\n\n")
  } else if (raw_or_LFQ == "raw") {
    cat("\n -- raw data are used --\n\n")
  }

  if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
  if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
  if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
  if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
  if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}
  if (name_column_groups == "thesearethesamplenamesused") {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_groups, thanks!")}

  if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
    GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
  }
  if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) != 2) {stop("To perform a Fold Change analyses, you must have exactly 2 groups!")}

  if (is.null(control_group)) {
    cat(paste0('\n"', levels(pull(GCPlist$sampleINFO, name_column_groups))[1], '" has been used as control group. If that is not fine for you, specify which group you want in the argument control_group\n\n'))
    second_to_first_ratio <- TRUE
  } else {
    if (length(control_group)!=1) {stop("control_group must be a character of length 1")}
    if (!is.character(control_group)) {stop("control_group must be a character of length 1")}
    if (is.na(control_group)) {stop("control_group must be a character of length 1, not a NA")}
    if (!control_group%in%levels(pull(GCPlist$sampleINFO, name_column_groups))) {stop("control_group must be a group contained in the column indicated by name_column_groups of the sampleINFO table")}
    if (control_group==levels(pull(GCPlist$sampleINFO, name_column_groups))[1]) {
      second_to_first_ratio <- TRUE
    } else {
      second_to_first_ratio <- FALSE
    }
  }

  if (length(paired)!=1) {stop("paired must be exclusively TRUE or FALSE")}
  if (!is.logical(paired)) {stop("paired must be exclusively TRUE or FALSE")}
  if (is.na(paired)) {stop("paired must be exclusively TRUE or FALSE")}

  if (paired) {
    length_vect_fact1 <- length(which(pull(GCPlist$sampleINFO, name_column_groups) == levels(pull(GCPlist$sampleINFO, name_column_groups))[1]))
    length_vect_fact2 <- length(which(pull(GCPlist$sampleINFO, name_column_groups) == levels(pull(GCPlist$sampleINFO, name_column_groups))[2]))

    if (length_vect_fact1 != length_vect_fact2) {stop("You cannot perform a paired Fold Change on groups that don't have the same number of observations")}
  }

  if (length(are_log_transf)!=1) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (!is.logical(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (is.na(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}

  if (length(log_base)!=1) {stop("log_base must be a number")}
  if (is.na(log_base)) {stop("log_base must be a number")}
  if (!is.numeric(log_base)) {stop("log_base must be a number")}


  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw, name_first_column = "thesearethesamplenamesused")
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ, name_first_column = "thesearethesamplenamesused")
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }

  if (any(map_lgl(df_intensities, ~ any(is.na(.))))) {stop("there are some missing values in the data")}


  df_intensities_wg <- add_column(df_intensities,
                                  allwiththis = factor(NA, levels = levels(pull(GCPlist$sampleINFO, name_column_groups))),
                                  .after = 1)
  colnames(df_intensities_wg)[2] <- name_column_groups

  for (i in 1:length(pull(df_intensities_wg, 1))) {
    df_intensities_wg[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
  }


  the_FC_table <- gentab_FC(df = df_intensities_wg,
                            v = colnames(df_intensities)[-1],
                            f = name_column_groups,
                            second_to_first_ratio = second_to_first_ratio,
                            paired = paired,
                            are_log_transf = are_log_transf,
                            log_base = log_base,
                            filter_sign = FALSE)
  colnames(the_FC_table)[1] <- "protid"

  the_FC_table <- mutate(the_FC_table, FCcomparison = rep(ifelse(second_to_first_ratio,
                                                                 paste0(levels(pull(df_intensities_wg, name_column_groups))[2], " vs ", levels(pull(df_intensities_wg, name_column_groups))[1]),
                                                                 paste0(levels(pull(df_intensities_wg, name_column_groups))[1], " vs ", levels(pull(df_intensities_wg, name_column_groups))[2])),
                                                          nrow(the_FC_table)))

  colnames_noprotid <- colnames(the_FC_table)[which(colnames(the_FC_table)!="protid")]
  colnames_noprotid_present <- colnames_noprotid[which(colnames_noprotid %in% colnames(GCPlist$proteinINFO))]
  colnames_noprotid_notpresent <- colnames_noprotid[which(!colnames_noprotid %in% colnames(GCPlist$proteinINFO))]

  if (length(colnames_noprotid_present)>0) {
    cat("\nThe following columns were already present in the proteinINFO data frame. Please note that they have now been replaced!\n ")
    cat(paste0(colnames_noprotid_present, collapse = "\n "))
    cat("\n")

    GCPlist$proteinINFO <- GCPlist$proteinINFO[, which(!colnames(GCPlist$proteinINFO) %in% colnames_noprotid_present)]
  }

  if (length(colnames_noprotid_notpresent)>0) {
    cat("\nThe following columns have been added to the proteinINFO table:\n ")
    cat(paste0(colnames_noprotid_notpresent, collapse = "\n "))
    cat("\n")
  }

  cat("\n")


  GCPoutput <- GCPlist

  GCPoutput$proteinINFO <- left_join(x = GCPlist$proteinINFO, y = the_FC_table, by = "protid")


  return(GCPoutput)
}
