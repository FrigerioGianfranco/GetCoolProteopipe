#' Get Principal Component analysis.
#'
#' It performs a principal component analysis.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The principal component analysis will be performed only in the specified data.
#'
#' @return The GCPlist will be returned with the scores in the sampleINFO and the loadings in the proteinINFO.
#'
#' @export
GCP_PCA <- function(GCPlist, raw_or_LFQ = c("raw", "LFQ")) {

  checkGCPlist(GCPlist)

  if (!identical(tolower(raw_or_LFQ), c("raw", "lfq"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("raw", "lfq"))


  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw, name_first_column = "thesearethesamplenamesused")
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ, name_first_column = "thesearethesamplenamesused")
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }


  PCA_list <- GetFeatistics::getPCA(df = df_intensities,
                                    v = colnames(df_intensities)[-1])

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
