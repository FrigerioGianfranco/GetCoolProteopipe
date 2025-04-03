#' Exporting the GCPlist
#'
#' It exports a GPClist in a single table file in the current working directory.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param filename character. The name of the file to create. It must end with ".txt", ".csv", or ".xlsx"; and the file will be accordingly created of that format.
#' @param exportype exportype one of the following: "raw", "LFQ", "proteinINFO", "sampleINFO", "all". If "raw", a table with the raw intensities will be exported; if "LFQ", a table with the LFQ intensities will be exported; if "proteinINFO", a table with all the proteinINFO; if "sampleINFO", a table with all the sampleINFO; if "all", a table with combined quant_raw, quant_LFQ, and proteinINFO will be exported. Please note that the latter might be too big to be suitably exported as a '.xlsx' file: if so, export it as '.txt' or '.csv'.
#' @param intensity_indication logical. if TRUE and if exportype is either "raw" or "LFQ" the column names of samples will start with "Intensity " for raw intensities or with "LFQ Intensity " for LFQ intensities. Please note that this will happen anyway if exportype is "all".
#' @param protgenenames logical. If TRUE and if exportype is either "raw", "LFQ", or "proteinINFO", it will export the columns 'Protein IDs' (or 'Accession', see below), 'Protein names', and 'Gene names'.
#' @param protIDclean logical. If TRUE (and if protgenenames is TRUE, and if exportype is either "raw", "LFQ", or "proteinINFO") it exports the column 'Accession', which contains only the first code of each protein of the Protein IDs. If FALSE it exports the complete 'Protein IDs' column.
#' @param specific_columns character. If exportype is "raw", "LFQ", or "proteinINFO", you can specify here some additional columns to export from the proteinINFO table. Moreover, you can simply pass here a single element to export all the related result columns, in particular: "shap_test", "PC", "ttest", "FC", "ANOVA"; and also "all_stat" for all of those.
#'
#' @return Export the file in the current working directory.
#'
#' @importFrom writexl write_xlsx
#'
#' @export
ExportGCPlist <- function(GCPlist, filename = "GCP.xlsx", exportype = c("raw", "LFQ", "proteinINFO", "sampleINFO", "all"), intensity_indication = TRUE, protgenenames = TRUE, protIDclean = TRUE, specific_columns = NULL) {

  checkGCPlist(GCPlist)

  if (length(filename) != 1) {stop('filename must be a character of length 1')}
  if (is.na(filename)) {stop('filename must be a character of length 1, not a missing value!')}
  if (!is.character(filename)) {stop('filename must be a character')}
  if (!(endsWith(toupper(filename), ".TXT") | endsWith(toupper(filename), ".CSV") | endsWith(toupper(filename), ".XLSX"))) {stop('filename must end with ".txt", ".csv", or ".xlsx"')}

  if (!identical(toupper(exportype), c("RAW", "LFQ", "PROTEININFO", "SAMPLEINFO", "ALL"))) {
    if (length(exportype) != 1) {stop('exportype must be one of "raw", "LFQ", "proteinINFO", "sampleINFO", "all"')}
    if (is.na(exportype)) {stop('exportype must be one of "raw", "LFQ", "proteinINFO", "sampleINFO", "all"')}
  }
  exportype <- toupper(exportype)
  exportype <- match.arg(exportype, c("RAW", "LFQ", "PROTEININFO", "SAMPLEINFO", "ALL"))



  if (length(intensity_indication) != 1) {stop("intensity_indication must be either TRUE or FALSE")}
  if (!is.logical(intensity_indication)) {stop("intensity_indication must be either TRUE or FALSE")}
  if (is.na(intensity_indication)) {stop("intensity_indication must be either TRUE or FALSE")}

  if (length(protgenenames) != 1) {stop("protgenenames must be either TRUE or FALSE")}
  if (!is.logical(protgenenames)) {stop("protgenenames must be either TRUE or FALSE")}
  if (is.na(protgenenames)) {stop("protgenenames must be either TRUE or FALSE")}

  if (length(protIDclean) != 1) {stop("protIDclean must be either TRUE or FALSE")}
  if (!is.logical(protIDclean)) {stop("protIDclean must be either TRUE or FALSE")}
  if (is.na(protIDclean)) {stop("protIDclean must be either TRUE or FALSE")}

  if (exportype %in% c(c("RAW", "LFQ", "PROTEININFO"))) {
    if (!is.null(specific_columns)) {
      if (!is.character(specific_columns)) {stop("specific_columns must be a character vector")}
      if (length(specific_columns) < 1) {
        specific_columns <- NULL
      } else {

        specific_columns_verified <- character()

        for (sc in specific_columns) {
          if (is.na(sc)) {
            stop("specific_columns must not contain NA")
          } else if (tolower(sc) == "all_stat") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(grepl("^shap_test_", colnames(GCPlist$proteinINFO)) |
                                                                                 grepl("^PC", colnames(GCPlist$proteinINFO)) |
                                                                                 grepl("^ttest_", colnames(GCPlist$proteinINFO)) |
                                                                                 colnames(GCPlist$proteinINFO)=="FC" | colnames(GCPlist$proteinINFO)=="logFC" | colnames(GCPlist$proteinINFO)=="FCcomparison" |
                                                                                 grepl("^ANOVA_", colnames(GCPlist$proteinINFO)))])
          } else if (tolower(sc) == "shap_test") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(grepl("^shap_test_", colnames(GCPlist$proteinINFO)))])
          } else if (tolower(sc) == "pc") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(grepl("^PC", colnames(GCPlist$proteinINFO)))])
          } else if (tolower(sc) == "ttest") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(grepl("^ttest_", colnames(GCPlist$proteinINFO)))])
          } else if (tolower(sc) == "fc") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(colnames(GCPlist$proteinINFO)=="FC" | colnames(GCPlist$proteinINFO)=="logFC" | colnames(GCPlist$proteinINFO)=="FCcomparison")])
          } else if (tolower(sc) == "anova") {
            specific_columns_verified <- c(specific_columns_verified,
                                           colnames(GCPlist$proteinINFO)[which(grepl("^ANOVA_", colnames(GCPlist$proteinINFO)))])
          } else if (!sc %in% colnames(GCPlist$proteinINFO)) {
            stop(paste0('the specific_column "', sc, '" does not exist among the columns of proteinINFO'))
          } else {
            specific_columns_verified <- c(specific_columns_verified, sc)
          }
        }

        specific_columns_verified <- unique(specific_columns_verified)
      }
    }
  }


  if (((toupper(exportype) == "RAW" | toupper(exportype) == "LFQ") & intensity_indication) | toupper(exportype) == "ALL") {
    colnames(GCPlist$quant_raw)[-which(colnames(GCPlist$quant_raw) == "protid")] <- paste0("Intensity ", colnames(GCPlist$quant_raw)[-which(colnames(GCPlist$quant_raw) == "protid")])
    colnames(GCPlist$quant_LFQ)[-which(colnames(GCPlist$quant_LFQ) == "protid")] <- paste0("LFQ Intensity ", colnames(GCPlist$quant_LFQ)[-which(colnames(GCPlist$quant_LFQ) == "protid")])
  }


  if (toupper(exportype) == "SAMPLEINFO") {
    table_to_export <- GCPlist$sampleINFO

  } else if (toupper(exportype) == "ALL") {

    table_to_export_1 <- left_join(x = GCPlist$proteinINFO, y = GCPlist$quant_raw, by = "protid", suffix = c("_INFO", "_raw"))
    table_to_export <- left_join(x = table_to_export_1, y = GCPlist$quant_LFQ, by = "protid", suffix = c("_INFO", "_LFQ"))


  } else {

    if (protgenenames == FALSE & is.null(specific_columns) & toupper(exportype) == "RAW") {
      table_to_export <- GCPlist$quant_raw
    } else if (protgenenames == FALSE & is.null(specific_columns) &  toupper(exportype) == "LFQ") {
      table_to_export <- GCPlist$quant_LFQ
    } else if (protgenenames == FALSE & is.null(specific_columns) &  toupper(exportype) == "PROTEININFO") {
      table_to_export <- GCPlist$proteinINFO
    } else {

      colproteinINFO_to_export <- c("protid")

      if (protgenenames) {

        if (protIDclean) {
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Accession")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Accession' column")}
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Protein names")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Protein names' column")}
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Gene names")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Gene names' column")}

          colproteinINFO_to_export <- c(colproteinINFO_to_export, "Accession", "Protein names", "Gene names")
        } else {
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Protein IDs")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Protein IDs' column")}
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Protein names")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Protein names' column")}
          if (length(which(colnames(GCPlist[["proteinINFO"]]) == "Gene names")) != 1) {stop("The proteinINFO dataframe of the GCPlist must have exactly one 'Gene names' column")}

          colproteinINFO_to_export <- c(colproteinINFO_to_export, "Protein IDs", "Protein names", "Gene names")
        }
      }


      if (!is.null(specific_columns)) {
        colproteinINFO_to_export <- c(colproteinINFO_to_export, specific_columns_verified)
      }

      colproteinINFO_to_export <- unique(colproteinINFO_to_export)


      if (toupper(exportype) == "RAW") {

        proteinINFO_fil <- GCPlist$proteinINFO[which(GCPlist$proteinINFO$protid%in%GCPlist$quant_raw$protid),]

        table_to_export <- left_join(x = select(proteinINFO_fil, all_of(colproteinINFO_to_export)), y = GCPlist$quant_raw, by = "protid", suffix = c("_INFO", "_raw"))

      } else if (toupper(exportype) == "LFQ") {

        proteinINFO_fil <- GCPlist$proteinINFO[which(GCPlist$proteinINFO$protid%in%GCPlist$quant_LFQ$protid),]

        table_to_export <- left_join(x = select(proteinINFO_fil, all_of(colproteinINFO_to_export)), y = GCPlist$quant_LFQ, by = "protid", suffix = c("_INFO", "_LFQ"))

      } else if (toupper(exportype) == "PROTEININFO") {

        table_to_export <- select(GCPlist$proteinINFO, all_of(colproteinINFO_to_export))

      }
    }
  }



  if (endsWith(toupper(filename), ".TXT")) {
    write_tsv(table_to_export, filename)
  } else if (endsWith(toupper(filename), ".CSV")) {
    write_csv(table_to_export, filename)
  } else if (endsWith(toupper(filename), ".XLSX")) {
    write_xlsx(table_to_export, filename)
  }
}
