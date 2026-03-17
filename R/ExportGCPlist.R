#' Exporting the GCPlist
#'
#' It exports a GPClist in a single table file in the current working directory.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function or a list of such lists. If the second is passed, the following filename argument must end with ".xlsx" as a sheet for each of those list will be created.
#' @param filename character. The name of the file to create. It must end with ".txt", ".csv", or ".xlsx"; and the file will be accordingly created of that format.
#' @param exportype exportype one of the following: "intensities", "proteinINFO", "sampleINFO", "all", "sheets". If "intensities", a table with the intensities will be exported; if "proteinINFO", a table with all the proteinINFO; if "sampleINFO", a table with all the sampleINFO; if "all", a table with combined intensities and proteinINFO will be exported; if "sheets", an Excel table with 3 sheets will be exported. Please note that "all" might be too big to be suitably exported as a '.xlsx' file: if so, export it as '.txt' or '.csv'.
#' @param intensity_indication logical. if TRUE and if exportype is "intensities", the column names of samples will start with "Intensity ". Please note that this will happen anyway if exportype is "all".
#' @param protgenenames logical. If TRUE and if exportype is "intensities" or "proteinINFO", it will export the columns 'Protein IDs' (or 'Accession', see below), 'Protein names', and 'Gene names'.
#' @param protIDclean logical. If TRUE (and if protgenenames is TRUE, and if exportype is "intensities" or "proteinINFO"), it exports the column 'Accession', which contains only the first code of each protein of the Protein IDs. If FALSE it exports the complete 'Protein IDs' column.
#' @param specific_columns character. If exportype is "intensities" or "proteinINFO", you can specify here some additional columns to export from the proteinINFO table. Moreover, you can simply pass here a single element to export all the related result columns, in particular: "shap_test", "PC", "ttest", "FC", "ANOVA", "presence_group"; and also "all_stat" for all of those.
#'
#' @return Export the file in the current working directory.
#'
#' @examples
#' \dontrun{
#'
#' # exporting a single GCPlist:
#'
#' ExportGCPlist(GCPlist14, "GCPlist14_single_export.txt")
#'
#'
#' # exporting multiple GCPlists:
#'
#' GCPlist_of_lists <- list(asimported = GCPlist00,
#'                          `before processing` = GCPlist05,
#'                          `post processing` = GCPlist11,
#'                          with_statistics = GCPlist14,
#'                          with_mock_ANOVA = GCPlist14a2,
#'                          with_mock_batchcorr = GCPlist14b2)
#'
#' ExportGCPlist(GCPlist_of_lists, "GCPlists_multiple_export.xlsx")
#'
#'
#' }
#'
#'
#' @importFrom writexl write_xlsx
#'
#' @export
ExportGCPlist <- function(GCPlist, filename = "GCP.xlsx", exportype = c("intensities", "proteinINFO", "sampleINFO", "all", "sheets"), intensity_indication = TRUE, protgenenames = TRUE, protIDclean = TRUE, specific_columns = "all_stat") {

  if (!is.list(GCPlist)) {stop("GCPlist must be a list")}
  if (length(GCPlist)<1) {stop("GCPlist is empty!")}
  if (is.data.frame(GCPlist[[1]])) {
    checkGCPlist(GCPlist)

    list_of_lists <- FALSE
  } else {

    if (!all(map_lgl(GCPlist, is.list))) {stop("GCPlist must be either a single GCPlist, or a list of such GCPlists")}

    for (ii in 1:length(GCPlist)) {
      checkGCPlist(GCPlist[[ii]])
      if (is.null(names(GCPlist)[ii])) {
        names(GCPlist)[ii] <- paste0("sheet", ifelse(ii<10, paste0("0", as.character(ii)), as.character(ii)))
      } else if (is.na(names(GCPlist)[ii])) {
        names(GCPlist)[ii] <- paste0("sheet", ifelse(ii<10, paste0("0", as.character(ii)), as.character(ii)))
      } else if (names(GCPlist)[ii] == "") {
        names(GCPlist)[ii] <- paste0("sheet", ifelse(ii<10, paste0("0", as.character(ii)), as.character(ii)))
      }
    }

    list_of_lists <- TRUE
  }

  if (length(filename) != 1) {stop('filename must be a character of length 1')}
  if (is.na(filename)) {stop('filename must be a character of length 1, not a missing value!')}
  if (!is.character(filename)) {stop('filename must be a character')}
  if (!(endsWith(toupper(filename), ".TXT") | endsWith(toupper(filename), ".CSV") | endsWith(toupper(filename), ".XLSX"))) {stop('filename must end with ".txt", ".csv", or ".xlsx"')}

  if (list_of_lists) {
    if (!endsWith(toupper(filename), ".XLSX")) {stop('Since you have a list of lists, filename must end with ".xlsx"')}
  }

  if (!identical(toupper(exportype), c("INTENSITIES", "PROTEININFO", "SAMPLEINFO", "ALL", "SHEETS"))) {
    if (length(exportype) != 1) {stop('exportype must be one of "intensities", "proteinINFO", "sampleINFO", "all", "sheets"')}
    if (is.na(exportype)) {stop('exportype must be one of "intensities", "proteinINFO", "sampleINFO", "all", "sheets"')}
  }
  exportype <- toupper(exportype)
  exportype <- match.arg(exportype, c("INTENSITIES", "PROTEININFO", "SAMPLEINFO", "ALL", "SHEETS"))

  if (toupper(exportype) == "SHEETS" & !endsWith(toupper(filename), ".XLSX")) {stop('If you select exportype as "sheets", the filename must end in ".xlsx"')}
  if (list_of_lists) {
    if (toupper(exportype) == "SHEETS") {stop('Since you have a list of lists, exportype cannot be "sheets"')}
  }

  if (length(intensity_indication) != 1) {stop("intensity_indication must be either TRUE or FALSE")}
  if (!is.logical(intensity_indication)) {stop("intensity_indication must be either TRUE or FALSE")}
  if (is.na(intensity_indication)) {stop("intensity_indication must be either TRUE or FALSE")}

  if (length(protgenenames) != 1) {stop("protgenenames must be either TRUE or FALSE")}
  if (!is.logical(protgenenames)) {stop("protgenenames must be either TRUE or FALSE")}
  if (is.na(protgenenames)) {stop("protgenenames must be either TRUE or FALSE")}

  if (length(protIDclean) != 1) {stop("protIDclean must be either TRUE or FALSE")}
  if (!is.logical(protIDclean)) {stop("protIDclean must be either TRUE or FALSE")}
  if (is.na(protIDclean)) {stop("protIDclean must be either TRUE or FALSE")}


  if (toupper(exportype) == "INTENSITIES") {
    cat("\n -- the intensities are being exported")

    if (protgenenames & protIDclean) {
      cat(",\n    preceded by 'Accession', 'Protein names', and 'Gene names'")
    } else if (protgenenames & !protIDclean) {
      cat(",\n    preceded by 'Protein IDs', 'Protein names', and 'Gene names'")
    }

    if (!is.null(specific_columns)) {
      cat(",\n    and followed by columns of results if present")
    }

    cat(" --\n\n")

  } else if (toupper(exportype) == "PROTEININFO") {
    cat("\n -- proteinINFO are being exported --\n\n")
  } else if (toupper(exportype) == "SAMPLEINFO") {
    cat("\n -- sampleINFO are being exported --\n\n")
  } else if (toupper(exportype) == "ALL") {
    cat("\n -- all of intensities and proteinINFO are being exported --\n\n")
  } else if (toupper(exportype) == "SHEETS") {
    cat("\n -- The tables intensities, proteinINFO, and sampleINFO of this list are being exported --\n\n")
  }



  things_to_do_to_each_GCPlist <- function(GCPlist, exportype, intensity_indication, protgenenames, protIDclean, specific_columns) {
    if (exportype %in% c("INTENSITIES", "PROTEININFO", "SHEETS")) {
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
                                                                                   endsWith(colnames(GCPlist$proteinINFO), "_FC") | endsWith(colnames(GCPlist$proteinINFO), "_logFC") |
                                                                                   grepl("^ANOVA_", colnames(GCPlist$proteinINFO)) |
                                                                                   grepl("^present_", colnames(GCPlist$proteinINFO)) | grepl("^combination_", colnames(GCPlist$proteinINFO)))])
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
                                             colnames(GCPlist$proteinINFO)[which(colnames(GCPlist$proteinINFO)=="FC" | colnames(GCPlist$proteinINFO)=="logFC" | colnames(GCPlist$proteinINFO)=="FCcomparison" |
                                                                                   endsWith(colnames(GCPlist$proteinINFO), "_FC") | endsWith(colnames(GCPlist$proteinINFO), "_logFC"))])
            } else if (tolower(sc) == "anova") {
              specific_columns_verified <- c(specific_columns_verified,
                                             colnames(GCPlist$proteinINFO)[which(grepl("^ANOVA_", colnames(GCPlist$proteinINFO)))])
            } else if (tolower(sc) == "presence_group") {
              specific_columns_verified <- c(specific_columns_verified,
                                             colnames(GCPlist$proteinINFO)[which(colnames(GCPlist$proteinINFO)=="present_" | colnames(GCPlist$proteinINFO)=="combination_")])
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


    if ((toupper(exportype) == "INTENSITIES" & intensity_indication) | toupper(exportype) == "ALL"| toupper(exportype) == "SHEETS") {
      colnames(GCPlist$intensities)[-which(colnames(GCPlist$intensities) == "protid")] <- paste0("Intensity ", colnames(GCPlist$intensities)[-which(colnames(GCPlist$intensities) == "protid")])
    }

    if (toupper(exportype) == "SAMPLEINFO") {
      table_to_export <- GCPlist$sampleINFO

    } else if (toupper(exportype) == "ALL") {

      table_to_export <- left_join(x = GCPlist$proteinINFO, y = GCPlist$intensities, by = "protid", suffix = c("_INFO", "_intensities"))

    } else {

      if (protgenenames == FALSE & is.null(specific_columns) & toupper(exportype) == "INTENSITIES") {
        table_to_export <- GCPlist$intensities
      } else if (protgenenames == FALSE & is.null(specific_columns) & toupper(exportype) == "PROTEININFO") {
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


        if (toupper(exportype) == "INTENSITIES") {

          proteinINFO_fil <- GCPlist$proteinINFO[which(GCPlist$proteinINFO$protid%in%GCPlist$intensities$protid),]

          table_to_export <- left_join(x = select(proteinINFO_fil, all_of(colproteinINFO_to_export)), y = GCPlist$intensities, by = "protid", suffix = c("_INFO", "_intensities"))

          if (!is.null(specific_columns)) {
            table_to_export <- relocate(table_to_export, all_of(specific_columns_verified), .after = last_col())
          }

        } else if (toupper(exportype) == "PROTEININFO") {

          table_to_export <- select(GCPlist$proteinINFO, all_of(colproteinINFO_to_export))

        }
      }
    }

    if (toupper(exportype) != "SHEETS") {
      return(table_to_export)
    } else {
      GCPout <- GCPlist

      GCPout$proteinINFO <- select(GCPout$proteinINFO, all_of(colproteinINFO_to_export))

      return(GCPout)
    }
  }


  if (!list_of_lists) {

    what_to_export <- things_to_do_to_each_GCPlist(GCPlist = GCPlist, exportype = exportype, intensity_indication = intensity_indication, protgenenames = protgenenames, protIDclean = protIDclean, specific_columns = specific_columns)

  } else {

    what_to_export <- vector(mode = "list", length = length(GCPlist))
    names(what_to_export) <- names(GCPlist)

    for (ll in 1:length(GCPlist)) {
      what_to_export[[ll]] <- things_to_do_to_each_GCPlist(GCPlist = GCPlist[[ll]], exportype = exportype, intensity_indication = intensity_indication, protgenenames = protgenenames, protIDclean = protIDclean, specific_columns = specific_columns)
    }

  }


  if (endsWith(toupper(filename), ".TXT")) {
    write_tsv(what_to_export, filename)
  } else if (endsWith(toupper(filename), ".CSV")) {
    write_csv(what_to_export, filename)
  } else if (endsWith(toupper(filename), ".XLSX")) {
    write_xlsx(what_to_export, filename)
  }
}
