#' Testing the normality
#'
#' It performs the Shapiro-Wilk test and print out a density plot and a qqplot.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param print_DensityPlot logical. If TRUE, a PDF with the density plots will be created in the current working directory.
#' @param print_QQPlot logical. If TRUE, a PDF with the Q-Q plots will be created in the current working directory.
#' @param print_only_these_protid NULL or a character vector of protid. Only the density/q-q-plot of the specified protid will be exported in the pdf files.
#' @param print_only_the_first NULL or a numeric integer. Only the density/q-q-plot of the those first proteins will be exported in the pdf files (suggested as usually with more than a thousand proteins, the file will be too big).
#'
#' @return The GCPlist with the results of the Shapiro-Wilk test in the proteinINFO table.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist05 <- GCP_TestNormality(GCPlist04)
#'
#' }
#'
#'
#' @export
GCP_TestNormality <- function(GCPlist, print_DensityPlot = FALSE, print_QQPlot = FALSE, print_only_these_protid = NULL, print_only_the_first = 10) {

  checkGCPlist(GCPlist)

  if (!is.logical(print_DensityPlot)) {stop("print_DensityPlot must be either TRUE or FALSE")}
  if (length(print_DensityPlot) != 1) {stop("print_DensityPlot must be either TRUE or FALSE")}
  if (is.na(print_DensityPlot)) {stop("print_DensityPlot must be either TRUE or FALSE")}

  if (!is.logical(print_QQPlot)) {stop("print_QQPlot must be either TRUE or FALSE")}
  if (length(print_QQPlot) != 1) {stop("print_QQPlot must be either TRUE or FALSE")}
  if (is.na(print_QQPlot)) {stop("print_QQPlot must be either TRUE or FALSE")}

  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities)

  if (!is.null(print_only_these_protid)) {
    if (!is.character(print_only_these_protid)) {stop("if not NULL, print_only_these_protid must be a character vector")}
    if (any(is.na(print_only_these_protid))) {stop("if not NULL, print_only_these_protid must not contain NAs")}
    if (!all(print_only_these_protid %in% colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])) {stop("if not NULL, print_only_these_protid must contain names of columns of the data intensities")}
  }

  if (!is.null(print_only_the_first)) {
    if (length(print_only_the_first)!=1) {stop("if not NULL, print_only_the_first must be a numeric of lenght 1")}
    if (!is.numeric(print_only_the_first)) {stop("if not NULL, print_only_the_first must be a numeric of lenght 1")}
    if (print_only_the_first<1) {stop("if not NULL, print_only_the_first must be a numeric integer greater or equal to 1")}
    if (print_only_the_first != as.integer(print_only_the_first)) {stop("if not NULL, print_only_the_first must be a numeric integer greater or equal to 1")}
    if (is.na(print_only_the_first)) {stop("if not NULL, print_only_the_first must be a numeric of lenght 1, not a missing value")}
  }

  if (!is.null(print_only_these_protid)) {
    prot_to_use <- print_only_these_protid
  } else if (!is.null(print_only_the_first)) {
    if (print_only_the_first > length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])) {
      print_only_the_first <- length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])
    }
    prot_to_use <- colnames(df_intensities)[-which(colnames(df_intensities) == "samples")][1:print_only_the_first]
  } else {
    prot_to_use <- colnames(df_intensities)[-which(colnames(df_intensities) == "samples")]
  }

  if (print_DensityPlot) {

    DensityPlots <- GetFeatistics::test_normality_density_plot(df = df_intensities, v = prot_to_use)

    GetFeatistics::export_figures(plots = DensityPlots,
                                  exprtname_figures = paste0("DensityPlots_", intensities),
                                  exprt_fig_type = "pdf")

  }

  if (print_QQPlot) {

    QQPlots <- GetFeatistics::test_normality_q_q_plot(df = df_intensities, v = prot_to_use)

    GetFeatistics::export_figures(plots = QQPlots,
                                  exprtname_figures = paste0("QQPlots_", intensities),
                                  exprt_fig_type = "pdf")
  }


  cat("\n_____\n")

  protid_with_at_least_3_notNA <- character()
  for (a in colnames(df_intensities)[-which(colnames(df_intensities) == "samples")]) {
    if (sum(!is.na(pull(df_intensities, a)))>=3) {
      protid_with_at_least_3_notNA <- c(protid_with_at_least_3_notNA, a)
    }
  }

  if (length(protid_with_at_least_3_notNA)!=length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])) {
    cat(paste0("There must be at least 3 non-missing values to perform a Shapiro-Wilk test:\nfor this reason, it is not being performed on ", length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])-length(protid_with_at_least_3_notNA), " proteins.\n\n"))
  }



  Shapiro_table <- GetFeatistics::test_normality_Shapiro_table(df = df_intensities, v = protid_with_at_least_3_notNA)


  cat(paste0("According to the Shapiro-Wilk test,\n",
             " ", sum(Shapiro_table$normally_distributed), " out of ", nrow(Shapiro_table), " proteins (", round(mean(Shapiro_table$normally_distributed)*100, digits = 1) , "%) are normally distributed\n",
             " ", sum(Shapiro_table$normally_distributed==FALSE), " out of ", nrow(Shapiro_table), " proteins (", round(mean(Shapiro_table$normally_distributed==FALSE)*100, digits = 1) , "%) are not normally distributed\n"))
  cat("_____\n")

  colnames(Shapiro_table)[1] <- "protid"
  colnames(Shapiro_table)[which(colnames(Shapiro_table) != "protid" & !grepl("^shap_test_", colnames(Shapiro_table)))] <- paste0("shap_test_", colnames(Shapiro_table)[which(colnames(Shapiro_table) != "protid" & !grepl("^shap_test_", colnames(Shapiro_table)))])


  colnames_noprotid <- colnames(Shapiro_table)[which(colnames(Shapiro_table)!="protid")]
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

  GCPoutput$proteinINFO <- left_join(x = GCPlist$proteinINFO, y = Shapiro_table, by = "protid")


  return(GCPoutput)

}
