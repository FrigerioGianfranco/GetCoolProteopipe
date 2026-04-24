#' performing a t-test on data
#'
#' It performs a t-test analyses on the proteins intensities.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups. Since this is a t-test, there must be exactly two groups.
#' @param paired logical. If FALSE it performs non-paired t-tests. If TRUE it performs paired t-tests.
#' @param FDR logical. If TRUE, after performing the t-tests, it also correct p-values across the different proteins with a false discovery rate multiple comparison correction (method "fdr" of the function p.adjust).
#' @param pcutoff a numeric of length 1, must be between 0 and 1. The difference between groups will be reported only if the p-values is below the cut-off reported here.
#'
#' @return The GCPlist with the results of the t-test added to the proteinINFO data frame.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist13 <- GCP_ttest(GCPlist = GCPlist12,
#'                        paired = TRUE,
#'                        FDR = TRUE,
#'                        pcutoff = 0.05)
#'
#' }
#'
#'
#'
#' @export
GCP_ttest <- function(GCPlist, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
                      paired = FALSE, FDR = TRUE, pcutoff = 0.05) {

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
  if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) != 2) {
    stop(paste0("To perform a t-test, you must have exactly 2 groups! The groups in '", name_column_groups, "' are:\n", paste0(levels(pull(GCPlist$sampleINFO, name_column_groups)), collapse = "\n")))
  }


  if (length(paired)!=1) {stop("paired must be exclusively TRUE or FALSE")}
  if (!is.logical(paired)) {stop("paired must be exclusively TRUE or FALSE")}
  if (is.na(paired)) {stop("paired must be exclusively TRUE or FALSE")}

  if (paired) {
    length_vect_fact1 <- length(which(pull(GCPlist$sampleINFO, name_column_groups) == levels(pull(GCPlist$sampleINFO, name_column_groups))[1]))
    length_vect_fact2 <- length(which(pull(GCPlist$sampleINFO, name_column_groups) == levels(pull(GCPlist$sampleINFO, name_column_groups))[2]))

    if (length_vect_fact1 != length_vect_fact2) {stop("You cannot perform a paired t-test on groups that don't have the same number of observations")}
  }

  if (length(FDR)!=1) {stop("FDR must be exclusively TRUE or FALSE")}
  if (!is.logical(FDR)) {stop("FDR must be exclusively TRUE or FALSE")}
  if (is.na(FDR)) {stop("FDR must be exclusively TRUE or FALSE")}

  if (length(pcutoff)!=1) {stop("pcutoff must be a single number between 0 and 1")}
  if (is.na(pcutoff)) {stop("pcutoff must be a single number between 0 and 1, and not a missing value")}
  if (!is.numeric(pcutoff)) {stop("pcutoff must be a single number between 0 and 1")}
  if (pcutoff>1 | pcutoff<0) {stop("pcutoff must be a single number between 0 and 1")}


  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities, name_first_column = "thesearethesamplenamesused")

  if (any(map_lgl(df_intensities, ~ any(is.na(.))))) {stop("there are some missing values in the data")}

  df_intensities_wg <- add_column(df_intensities,
                                  allwiththis = factor(NA, levels = levels(pull(GCPlist$sampleINFO, name_column_groups))),
                                  .after = 1)
  colnames(df_intensities_wg)[2] <- name_column_groups

  for (i in 1:length(pull(df_intensities_wg, 1))) {
    df_intensities_wg[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
  }


  the_ttest_table <- GetFeatistics::gentab_P.t.test(df = df_intensities_wg,
                                                    v = colnames(df_intensities)[-1],
                                                    f = name_column_groups,
                                                    paired = paired,
                                                    FDR = FDR,
                                                    cutPval = FALSE,
                                                    groupdiff = TRUE,
                                                    pcutoff = pcutoff,
                                                    filter_sign = FALSE)
  colnames(the_ttest_table)[1] <- "protid"
  colnames(the_ttest_table)[which(colnames(the_ttest_table)!="protid")] <- paste0("ttest_", colnames(the_ttest_table)[which(colnames(the_ttest_table)!="protid")])

  colnames_noprotid <- colnames(the_ttest_table)[which(colnames(the_ttest_table)!="protid")]
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

  GCPoutput$proteinINFO <- left_join(x = GCPlist$proteinINFO, y = the_ttest_table, by = "protid")


  return(GCPoutput)
}
