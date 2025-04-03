#' Filter proteins based on missing values per condition
#'
#' In the quant_raw and the quantLFQ data intensity, it keeps only proteins that are not missing values for a defined ratio in at least one of the specified conditions (sample groups).
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param ratio numeric, between 0 an 1 (included). Ratio of non-missing values, per sample group, wanted to pass this filtration step (the higher this ratio, the less proteins will be kept).
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#'
#' @return a GCPlist list with filtered tables.
#'
#' @export
GCP_FilterNAperCondition <- function(GCPlist, ratio = 0.5, name_column_groups = NULL) {

  checkGCPlist(GCPlist)

  if (length(ratio)!=1) {stop("ratio must be a numeric of length 1")}
  if (!is.numeric(ratio)) {stop("ratio must be a numeric of length 1")}
  if (is.na(ratio)) {stop("ratio must be a numeric of length 1, not a missing value")}
  if (ratio<0 | ratio>1) {stop("ratio must be between 0 and 1 (included)")}

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (any(pull(GCPlist$sampleINFO, name_column_groups)%in%c("prot", "so_to_keep"))) {stop("Please, don't call any group 'prot' or 'so_to_keep'")}
    if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}
  }



  if (is.null(name_column_groups)) {

    if ("allwiththis" %in% colnames(GCPlist$sampleINFO)) {stop("Please, don't call a column of the sampleINFO data frame 'allwiththis' as I need to create one with this name now")}

    cat("\n")
    cat("Since the name_column_groups was not specified, the filtration is performed taking into account all the samples")

    GCPlist$sampleINFO <- mutate(GCPlist$sampleINFO, allwiththis = as.factor("theOnlyGroup"))

    name_column_groups <- "allwiththis"
  }


  table_info_to_keep_raw <- tibble(prot = pull(GCPlist$quant_raw, "protid"))
  table_info_to_keep_raw[,unique(as.character(pull(GCPlist$sampleINFO, name_column_groups)))] <- NA
  table_info_to_keep_raw[,"so_to_keep"] <- NA

  table_info_to_keep_LFQ <- tibble(prot = pull(GCPlist$quant_LFQ, "protid"))
  table_info_to_keep_LFQ[,unique(as.character(pull(GCPlist$sampleINFO, name_column_groups)))] <- NA
  table_info_to_keep_LFQ[,"so_to_keep"] <- NA

  for (g in unique(as.character(pull(GCPlist$sampleINFO, name_column_groups)))) {

    name_samples <- pull(GCPlist$sampleINFO, 1)[which(pull(GCPlist$sampleINFO, name_column_groups) == g)]

    for (ir in 1:length(pull(GCPlist$quant_raw, 1))) {

      this_vector_raw <- as.vector(GCPlist$quant_raw[ir, name_samples], mode = "numeric")

      if (length(which(!is.na(this_vector_raw)))/length(this_vector_raw) >= ratio) {
        table_info_to_keep_raw[ir,g] <- TRUE
      } else {
        table_info_to_keep_raw[ir,g] <- FALSE
      }

      if (g == unique(as.character(pull(GCPlist$sampleINFO, name_column_groups)))[length(unique(as.character(pull(GCPlist$sampleINFO, name_column_groups))))]) {
        table_info_to_keep_raw[ir, "so_to_keep"] <- any(as.vector(table_info_to_keep_raw[ir, colnames(table_info_to_keep_raw)[which(!colnames(table_info_to_keep_raw)%in%c("prot", "so_to_keep"))]], mode = "logical"))
      }

    }

    for (il in 1:length(pull(GCPlist$quant_LFQ, 1))) {

      this_vector_LFQ <- as.vector(GCPlist$quant_LFQ[il, name_samples], mode = "numeric")

      if (length(which(!is.na(this_vector_LFQ)))/length(this_vector_LFQ) >= ratio) {
        table_info_to_keep_LFQ[il,g] <- TRUE
      } else {
        table_info_to_keep_LFQ[il,g] <- FALSE
      }

      if (g == unique(as.character(pull(GCPlist$sampleINFO, name_column_groups)))[length(unique(as.character(pull(GCPlist$sampleINFO, name_column_groups))))]) {
        table_info_to_keep_LFQ[il, "so_to_keep"] <- any(as.vector(table_info_to_keep_LFQ[il, colnames(table_info_to_keep_LFQ)[which(!colnames(table_info_to_keep_LFQ)%in%c("prot", "so_to_keep"))]], mode = "logical"))
      }

    }
  }


  GCPoutput <- GCPlist

  GCPoutput$quant_raw <- GCPoutput$quant_raw[which(pull(table_info_to_keep_raw, "so_to_keep")),]
  GCPoutput$quant_LFQ <- GCPoutput$quant_LFQ[which(pull(table_info_to_keep_LFQ, "so_to_keep")),]

  cat("\n")
  cat("_____________________\n")
  cat("For the quant_raw table:")

  if (name_column_groups != "allwiththis") {
    for (gr in colnames(table_info_to_keep_raw)[which(!colnames(table_info_to_keep_raw) %in% c("prot", "so_to_keep"))]) {
      cat("\n")
      cat(paste0("- ", length(which(pull(table_info_to_keep_raw, gr))), " out of ", length(pull(table_info_to_keep_raw)), " where suitable for ", gr, "."))
    }
  }
  cat("\n")
  cat(paste0("--> Overall, ", length(which(pull(table_info_to_keep_raw, "so_to_keep"))), " out of ", length(pull(table_info_to_keep_raw, "so_to_keep")), " have been kept in the quant_raw table."))
  cat("\n")
  cat("___________________\n")
  cat("\n")
  cat("For the quant_LFQ table:")
  if (name_column_groups != "allwiththis") {
    for (gl in colnames(table_info_to_keep_LFQ)[which(!colnames(table_info_to_keep_LFQ) %in% c("prot", "so_to_keep"))]) {
      cat("\n")
      cat(paste0("- ", length(which(pull(table_info_to_keep_LFQ, gl))), " out of ", length(pull(table_info_to_keep_LFQ)), " where suitable for ", gl, "."))
    }
  }
  cat("\n")
  cat(paste0("--> Overall, ", length(which(pull(table_info_to_keep_LFQ, "so_to_keep"))), " out of ", length(pull(table_info_to_keep_LFQ, "so_to_keep")), " have been kept in the quant_LFQ table."))
  cat("\n")
  cat("_____________________\n")


  return(GCPoutput)

}
