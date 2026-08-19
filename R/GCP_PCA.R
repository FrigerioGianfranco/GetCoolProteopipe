#' Get Principal Component analysis.
#'
#' It performs a principal component analysis.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#' @param center logical. Whether the variables should be shifted to be zero centered (as in the prcomp function).
#' @param scale. logical. whether the variables should be scaled to have unit variance before the analysis takes place (as in prcomp function).
#'
#' @return The GCPlist will be returned with the scores in the sampleINFO and the loadings in the proteinINFO.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist12 <- GCP_PCA(GCPlist = GCPlist11)
#'
#' }
#'
#'
#'
#' @export
GCP_PCA <- function(GCPlist, center = TRUE, scale. = FALSE) {

  checkGCPlist(GCPlist)

  if (length(center)!=1) {stop("center must be exclusively TRUE or FALSE")}
  if (!is.logical(center)) {stop("center must be exclusively TRUE or FALSE")}
  if (is.na(center)) {stop("center must be exclusively TRUE or FALSE")}

  if (length(scale.)!=1) {stop("scale. must be exclusively TRUE or FALSE")}
  if (!is.logical(scale.)) {stop("scale. must be exclusively TRUE or FALSE")}
  if (is.na(scale.)) {stop("scale. must be exclusively TRUE or FALSE")}


  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities, name_first_column = "thesearethesamplenamesused")

  PCA_list <- GetFeatistics::getPCA(df = df_intensities,
                                    v = colnames(df_intensities)[-1],
                                    center = center,
                                    scale. = scale.)

  GCPoutput <- GCPlist

  the_score_table <- PCA_list$df_with_scores_table[, colnames(PCA_list$df_with_scores_table)[which(!colnames(PCA_list$df_with_scores_table)%in%colnames(df_intensities)[-1])]]

  colnames(the_score_table)[1] <- colnames(GCPlist$sampleINFO)[1]
  colnames(the_score_table)[which(colnames(the_score_table)!=colnames(GCPlist$sampleINFO)[1])] <- paste0(colnames(the_score_table)[which(colnames(the_score_table)!=colnames(GCPlist$sampleINFO)[1])], "_scores")

  colnames_scoretab_noprotid <- colnames(the_score_table)[which(colnames(the_score_table)!=colnames(GCPlist$sampleINFO)[1])]
  colnames_scoretab_noprotid_present <- colnames_scoretab_noprotid[which(colnames_scoretab_noprotid %in% colnames(GCPlist$sampleINFO))]
  colnames_scoretab_noprotid_notpresent <- colnames_scoretab_noprotid[which(!colnames_scoretab_noprotid %in% colnames(GCPlist$sampleINFO))]

  if (length(colnames_scoretab_noprotid_present)>0) {
    cat("\nThe following columns were already present in the sampleINFO data frame. Please note that they have now been replaced!\n ")
    cat(paste0(colnames_scoretab_noprotid_present, collapse = "\n "))
    cat("\n")

    GCPlist$sampleINFO <- GCPlist$sampleINFO[, which(!colnames(GCPlist$sampleINFO) %in% colnames_scoretab_noprotid_present)]
  }

  if (length(colnames_scoretab_noprotid_notpresent)>0) {
    cat("\nThe following columns have been added to the sampleINFO table:\n ")
    cat(paste0(colnames_scoretab_noprotid_notpresent, collapse = "\n "))
    cat("\n")
  }


  GCPoutput$sampleINFO <- left_join(GCPlist$sampleINFO, the_score_table, by = colnames(GCPlist$sampleINFO)[1])


  the_loading_table <- PCA_list$dfv_with_loadings_table[, colnames(PCA_list$dfv_with_loadings_table)[which(!colnames(PCA_list$dfv_with_loadings_table)%in%colnames(GCPlist$proteinINFO)[-1])]]

  colnames(the_loading_table)[1] <- "protid"
  colnames(the_loading_table)[which(colnames(the_loading_table)!="protid")] <- paste0(colnames(the_loading_table)[which(colnames(the_loading_table)!="protid")], "_loadings")

  colnames_noprotid <- colnames(the_loading_table)[which(colnames(the_loading_table)!="protid")]
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

  GCPoutput$proteinINFO <- left_join(x = GCPlist$proteinINFO, y = the_loading_table, by = "protid")


  return(GCPoutput)
}
