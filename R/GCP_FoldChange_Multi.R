#' Performing a Fold Change analysis on data, for more than 2 groups
#'
#' It performs a Fold Change analyses on the proteins intensities, performing multiple pair comparisons. Please, be aware that the Fold Change analysis should be performed only on positive data! Indeed, all protein intensities should be positive, or unreliable results will be generated!
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function. IMPORTANT: All protein intensities should be positive, or unreliable results will be generated!
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups. 3 or more groups should be indicated here.
#' @param group_order NULL or a character. The ordered names of the groups, the first ones in order will be considered as denominator in the ratio of the Fold Change analyses. If NULL, the order of the levels of the factor will be considered.
#' @param paired logical. If FALSE it performs FC on mean of the two groups, for each pair. If TRUE it performs FC for each pair and then compute the mean.
#' @param are_log_transf logical. If the protein intensities are already log-transformed, specify here as TRUE, so the subtraction will be performed instead of the ratio.
#' @param log_base numeric of length 1. Specify here the base of the logarithm to calculate the logFC or, if are_log_transf is TRUE; the base of the logarithm that were used to transform the data.
#'
#' @return The GCPlist with the results of the Fold Change analyses added to the proteinINFO data frame.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist14a2 <- GCP_FoldChange_Multi(GCPlist = GCPlist14a1,
#'                                     name_column_groups = "Group_multi",
#'                                     paired = FALSE,
#'                                     are_log_transf = TRUE,
#'                                     log_base = 2)
#'
#' }
#'
#'
#'
#' @export
GCP_FoldChange_Multi <- function(GCPlist, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), group_order = NULL,
                                 paired = FALSE, are_log_transf = TRUE, log_base = getOption("GetCoolProteopipe.log_base")) {

  checkGCPlist(GCPlist)

  if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
  if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
  if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
  cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
  if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
  if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}
  if (name_column_groups == "thesearethesamplenamesused") {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_groups, thanks!")}

  if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
    GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
  }
  if (any(table(pull(GCPlist$sampleINFO, name_column_groups)) == 0)) {
    GCPlist$sampleINFO[,name_column_groups] <- droplevels(pull(GCPlist$sampleINFO, name_column_groups))
  }

  if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) <= 2) {
    stop(paste0("To perform these multiple Fold Change analyses, you should have more than 2 groups! The groups in '", name_column_groups, "' are:\n", paste0(levels(pull(GCPlist$sampleINFO, name_column_groups)), collapse = "\n")))
  }

  if (!is.null(group_order)) {
    if (length(group_order)!=length(levels(pull(GCPlist$sampleINFO, name_column_groups)))) {stop("the length of group_order must be the same of the groups")}
    if (any(is.na(group_order))) {stop("group_order must not contain NAs")}
    if (!is.character(group_order)) {stop("group_order must be a character")}
    if (any(duplicated(group_order))) {stop("group_order must not contain duplicated")}
    if (!all(group_order %in% levels(pull(GCPlist$sampleINFO, name_column_groups)))) {
      stop(paste0("group_order must contain the levels of the factor! The levels in '", name_column_groups, "' are:\n", paste0(levels(pull(GCPlist$sampleINFO, name_column_groups)), collapse = "\n")))
    }

    GCPlist$sampleINFO[,name_column_groups] <- factor(as.character(pull(GCPlist$sampleINFO, name_column_groups)), levels = group_order)
  }

  if (length(paired)!=1) {stop("paired must be exclusively TRUE or FALSE")}
  if (!is.logical(paired)) {stop("paired must be exclusively TRUE or FALSE")}
  if (is.na(paired)) {stop("paired must be exclusively TRUE or FALSE")}

  if (length(are_log_transf)!=1) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (!is.logical(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}
  if (is.na(are_log_transf)) {stop("are_log_transf must be exclusively TRUE or FALSE")}

  if (length(log_base)!=1) {stop("log_base must be a number")}
  if (is.na(log_base)) {stop("log_base must be a number")}
  if (!is.numeric(log_base)) {stop("log_base must be a number")}

  cat(paste0("\n -- The log_base considered is ", log_base, " --\n\n"))


  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities, name_first_column = "thesearethesamplenamesused")

  if (any(map_lgl(df_intensities, ~ any(is.na(.))))) {stop("there are some missing values in the data")}


  df_intensities_wg <- add_column(df_intensities,
                                  allwiththis = factor(NA, levels = levels(pull(GCPlist$sampleINFO, name_column_groups))),
                                  .after = 1)
  colnames(df_intensities_wg)[2] <- name_column_groups

  for (i in 1:length(pull(df_intensities_wg, 1))) {
    df_intensities_wg[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
  }


  the_FC_table_Multi <- gentab_FC_more_than2levels(df = df_intensities_wg,
                                                   v = colnames(df_intensities)[-1],
                                                   f = name_column_groups,
                                                   second_to_first_ratio = TRUE,
                                                   paired = paired,
                                                   are_log_transf = are_log_transf,
                                                   log_base = log_base,
                                                   only_on_positive = FALSE)

  colnames(the_FC_table_Multi)[1] <- "protid"




  colnames_noprotid <- colnames(the_FC_table_Multi)[which(colnames(the_FC_table_Multi)!="protid")]
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

  GCPoutput$proteinINFO <- left_join(x = GCPlist$proteinINFO, y = the_FC_table_Multi, by = "protid")


  return(GCPoutput)
}
