#' Know proteins shared by groups.
#'
#' It creates an additional column in the proteinINFO data frame containing where in which groups proteins are present.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups. The sample groups must be between at least 2.
#' @param raw_or_LFQ one of the following: "raw", "LFQ".
#'
#' @return The GCPlist with an additional column in the proteinINFO data frame.
#'
#' @export
GCP_ProteinsGrouped <- function(GCPlist, name_column_groups = NULL, raw_or_LFQ = getOption("GetCoolProteopipe.raw_or_LFQ")) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

    if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) < 2) {stop("The sample groups must be at least 2")}

  } else {
    stop("please, specifiy the name of the sample group column in the name_column_groups argoument")
  }

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


  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw, name_first_column = colnames(GCPlist$sampleINFO)[1])
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ, name_first_column = colnames(GCPlist$sampleINFO)[1])
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }

  df_intensities_w_groups <- left_join(x = GCPlist$sampleINFO[,unique(c(colnames(GCPlist$sampleINFO)[1], name_column_groups))], y = df_intensities, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_intensities"))

  list_for_ggven <- vector(mode = "list", length = length(levels(pull(GCPlist$sampleINFO, name_column_groups))))
  names(list_for_ggven) <- levels(pull(GCPlist$sampleINFO, name_column_groups))

  for (gr in levels(pull(GCPlist$sampleINFO, name_column_groups))) {

    df_intensities_w_groups_fill <- df_intensities_w_groups[which(pull(df_intensities_w_groups, name_column_groups) == gr), colnames(df_intensities_w_groups)[which(!colnames(df_intensities_w_groups)%in%c(colnames(GCPlist$sampleINFO)[1], name_column_groups))]]

    list_for_ggven[[gr]] <- character()

    for (a in colnames(df_intensities_w_groups_fill)) {
      if (any(!is.na(pull(df_intensities_w_groups_fill, a)))) {
        list_for_ggven[[gr]] <- c(list_for_ggven[[gr]], a)
      }
    }
  }



  GCPoutout <- GCPlist


  GCPoutout$proteinINFO[, paste0("present_", raw_or_LFQ, "_", levels(pull(GCPlist$sampleINFO, name_column_groups)))] <- as.logical(NA)
  GCPoutout$proteinINFO[, paste0("combination_", raw_or_LFQ, "_", name_column_groups)] <- as.character(NA)

  new_col_names <- c(paste0("present_", raw_or_LFQ, "_", levels(pull(GCPlist$sampleINFO, name_column_groups))), paste0("combination_", raw_or_LFQ, "_", name_column_groups))

  if (any(new_col_names %in% colnames(GCPlist$proteinINFO))) {
    warning(paste0('\nThe following columns were already present in proteinINFO, so now they have been replaced.\n "',
                   paste0(new_col_names[which(new_col_names %in% colnames(GCPlist$proteinINFO))], collapse = '", "'), '"'))
  }


  for (i in 1:nrow(GCPlist$proteinINFO)) {

    combining_presence <- character()

    for (gr in levels(pull(GCPlist$sampleINFO, name_column_groups))) {
      if (GCPoutout$proteinINFO$protid[i] %in% list_for_ggven[[gr]]) {
        GCPoutout$proteinINFO[i, paste0("present_", raw_or_LFQ, "_",  gr)] <- TRUE
        combining_presence <- c(combining_presence, gr)
      } else {
        GCPoutout$proteinINFO[i, paste0("present_", raw_or_LFQ, "_", gr)] <- FALSE
      }
    }

    GCPoutout$proteinINFO[i, paste0("combination_", raw_or_LFQ, "_",  name_column_groups)] <- paste0(combining_presence, collapse = "_")
  }


  all_combinations <- function(x) {
    combos <- unlist(
      lapply(1:length(x), function(i) {
        apply(combn(x, i), 2, function(y) paste(y, collapse = "_"))
      }),
      use.names = FALSE
    )
    c("", combos)
  }


  GCPoutout$proteinINFO[, paste0("combination_", raw_or_LFQ, "_",  name_column_groups)] <- factor(pull(GCPoutout$proteinINFO, paste0("combination_", raw_or_LFQ, "_",  name_column_groups)),
                                                                                                  levels = all_combinations(levels(pull(GCPlist$sampleINFO, name_column_groups))))




  return(GCPoutout)
}
