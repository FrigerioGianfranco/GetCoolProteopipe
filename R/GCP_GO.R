#' Performing Gene Ontology (GO) Enrichment Analyses
#'
#' It performs multiple Gene Ontology (GO) Enrichment Analyses for each pairwise comparison using the function enrichGO() from the clusterProfiler package. The genes used are the ones retrieved from the 'Gene names' column of proteinINFO.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups. If NULL, you will need to specify Pvalue_column and, possibly, FC_column, see below.
#' @param use_FDRPvalue logical. If name_column_groups is specified and if this is TRUE, this function will consider columns with FDR p-values in proteinINFO; otherwise, it will consider columns with the p-values in proteinINFO.
#' @param Pvalue_column NULL or character. If name_column_groups was not specified, you can directly pass here the column names of proteinINFO you want to consider as p-values.
#' @param Pvalue_cutoff numeric. P-values below this numner will be consider as statistically significant and relative gene names will be considered.
#' @param differential_analysis logical. If TRUE, it will perform two separate analyses for each pairwise condition considering the Fold Changes, see below;if FALSE; it will only perform an analyses for each pairwise condition.
#' @param use_logFC logical. If name_column_groups is specified and if this is TRUE, this function will consider columns with log FC in proteinINFO; otherwise, it will consider columns with FC in proteinINFO.
#' @param FC_cutoff numeric. significant FC higher than this value will be considered for the differential 'up' analysis, and significant FC lower than this value will be considered for the differential 'down' analysis.
#' @param background_from_data logical. If TRUE; the "universe" passed to enrichGO will be valid unique gene names from all the not-NA protein groups; if FALSE, all genes listed in the database will be used as background.
#' @param OrgDb One of "human", "mouse", or "rat"; or an OrgDb object to pass to clusterProfiler::enrichGO. If "human", the OrgDb "org.Hs.eg.db" will be used; if "mouse", "org.Mm.eg.db"; if "rat", "org.Rn.eg.db". In any case, you need to have that package installed; if you don't, please follow the instruction explained in the error message to do so.
#' @param keyType character of length 1. The keyType to pass to clusterProfiler::enrichGO. Considering that the values are taken from the 'Gene names' of proteinINFO, this should be "SYMBOL" as the defalut.
#' @param ont One of "ALL", "BP", "MF", or "CC". ont to pass to clusterProfiler::enrichGO.
#' @param pvalueCutoff numeric. pvalueCutoff to pass to clusterProfiler::enrichGO.
#' @param pAdjustMethod one of "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr". pAdjustMethod to pass to clusterProfiler::enrichGO.
#' @param qvalueCutoff numeric. qvalueCutoff to pass to clusterProfiler::enrichGO.
#' @param ... Further argument to pass to clusterProfiler::enrichGO.
#'
#' @return a list with three elements:
#'
#' - `considered_gene_names`: a list with all the gene names considered.
#'
#' - `SUMMARY_gene_names`: a tibble which is a summary of the numbers of gene names considered.
#'
#' - `GO_results`: A list with dataframes, each containing the results of the GO.
#'
#'
#' @examples
#' \dontrun{
#'
#'
#'
#' # specifying a name_column_groups:
#'
#' the_list_of_GO_results <- GCP_GO(GCPlist = GCPlist12,
#'                                  name_column_groups = "Condition",
#'                                  use_FDRPvalue = TRUE,
#'                                  Pvalue_cutoff = 0.05,
#'                                  differential_analysis = TRUE,
#'                                  use_logFC = TRUE,
#'                                  FC_cutoff = 0,
#'                                  background_from_data = TRUE,
#'                                  OrgDb = "mouse")
#'
#' # not specifying a name_column_groups:
#'
#' the_list_of_GO_results <- GCP_GO(GCPlist = GCPlist12,
#'                                  name_column_groups = NULL,
#'                                  Pvalue_column = c("ttest_V_vs_S_PvaluesFDR"),
#'                                  Pvalue_cutoff = 0.05,
#'                                  differential_analysis = TRUE,
#'                                  FC_column = c("V_vs_S_logFC"),
#'                                  FC_cutoff = 0,
#'                                  background_from_data = TRUE,
#'                                  OrgDb = "mouse")
#'
#'
#' }
#'
#' @importFrom clusterProfiler enrichGO
#'
#' @export
GCP_GO <- function(GCPlist,
                   name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
                   use_FDRPvalue = TRUE, Pvalue_column = NULL, Pvalue_cutoff = 0.05,
                   differential_analysis = TRUE, use_logFC = TRUE, FC_column = NULL, FC_cutoff = 0,
                   background_from_data = TRUE,
                   OrgDb = "human", keyType = "SYMBOL", ont = c("ALL", "BP", "MF", "CC"),
                   pvalueCutoff = 0.05, pAdjustMethod = "none", qvalueCutoff = 0.2, ...) {



  checkGCPlist(GCPlist)

  if (!"Gene names"%in% colnames(GCPlist$proteinINFO)) {stop("'Gene names' must be a column of the proteinINFO table of the GCPlist")}
  if (length(which(colnames(GCPlist$proteinINFO)=="Gene names"))!=1) {stop("'Gene names' must be unique column of the proteinINFO table of the GCPlist")}

  if (length(differential_analysis)!=1) {stop("differential_analysis must be TRUE or FALSE")}
  if (is.na(differential_analysis)) {stop("differential_analysis must be TRUE or FALSE")}
  if (!is.logical(differential_analysis)) {stop("differential_analysis must be TRUE or FALSE")}

  if (length(use_FDRPvalue)!=1) {stop("use_FDRPvalue must be TRUE or FALSE")}
  if (is.na(use_FDRPvalue)) {stop("use_FDRPvalue must be TRUE or FALSE")}
  if (!is.logical(use_FDRPvalue)) {stop("use_FDRPvalue must be TRUE or FALSE")}

  if (length(use_logFC)!=1) {stop("use_logFC must be TRUE or FALSE")}
  if (is.na(use_logFC)) {stop("use_logFC must be TRUE or FALSE")}
  if (!is.logical(use_logFC)) {stop("use_logFC must be TRUE or FALSE")}

  check_pattern_in_groups <- function(pattern, groups_considered) {
    if (any(grepl(pattern, groups_considered))) {
      stop(paste0("Unfortunately, the following groups cannot be considered as contain the pattern '", pattern , "' that cannot be considered:\n",
                  paste0(groups_considered[which(grepl(pattern, groups_considered))], collapse = "\n")))
    }
  }

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("if not NULL, name_column_groups must be a character of length 1")}
    if (!is.character(name_column_groups)) {stop("if not NULL, name_column_groups must be a character of length 1")}
    if (is.na(name_column_groups)) {stop("if not NULL, name_column_groups must be a character of length 1, not a missing value")}
    cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[, name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }
    groups_considered <- levels(pull(GCPlist$sampleINFO, name_column_groups))



    check_pattern_in_groups(pattern = "_vs_", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "ttest_", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "ANOVA_", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "_PvaluesFDR", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "_PvalueFDR", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "_Pvalues", groups_considered = groups_considered)
    check_pattern_in_groups(pattern = "_Pvalue", groups_considered = groups_considered)

    if (length(groups_considered) < 2) {
      stop(paste0("There are less than two groups in the column passed to name_column_groups '", name_column_groups, "'"))
    } else if (length(groups_considered) == 2) {
      statistica_test_to_consider <- "ttest"
      cat("\nSince there are two groups, the results of the ttest are considered\n")
    } else {
      statistica_test_to_consider <- "ANOVA"
      cat("\nSince there are more than two groups, the results of the ANOVA are considered\n")
    }

    all_possible_pairs <- as.vector(unlist(combn(groups_considered, 2, \(x) {c(paste(x[1], x[2], sep = "_vs_"), paste(x[2], x[1], sep = "_vs_"))})))

    if (!is.null(Pvalue_column)) {
      cat("\nWhat you passed to the argument Pvalue_column has been ignored as that information is taked from the name_column_groups you passed\n")
    }

    if (use_FDRPvalue) {
      all_possible_pairs_stat <- paste0(statistica_test_to_consider, "_", all_possible_pairs, ifelse(statistica_test_to_consider=="ttest", "_PvaluesFDR", "_PvalueFDR"))
    } else {
      all_possible_pairs_stat <- paste0(statistica_test_to_consider, "_", all_possible_pairs, ifelse(statistica_test_to_consider=="ttest", "_Pvalues", "_Pvalue"))
    }

    Pvalue_column <- all_possible_pairs_stat[which(all_possible_pairs_stat%in%colnames(GCPlist$proteinINFO))]

    if (length(Pvalue_column) < 1) {
      stop(paste0("Considering the name_column_groups '", name_column_groups, "' and that it contains ", length(groups_considered), " groups so ", statistica_test_to_consider, " is considered,\n",
                  "no suitable Pvalue columns are found in the proteinINFO table"))
    }

    if (differential_analysis) {

      check_pattern_in_groups(pattern = "_FC", groups_considered = groups_considered)
      check_pattern_in_groups(pattern = "_logFC", groups_considered = groups_considered)

      if (!is.null(FC_column)) {
        cat("\nWhat you passed to the argument FC_column has been ignored as that information is taked from the name_column_groups you passed\n")
      }

      if (use_logFC) {
        all_possible_pairs_FC <- paste0(all_possible_pairs, "_logFC")
      } else {
        all_possible_pairs_FC <- paste0(all_possible_pairs, "_FC")
      }

      FC_column <- all_possible_pairs_FC[which(all_possible_pairs_FC%in%colnames(GCPlist$proteinINFO))]

      if (length(FC_column) < 1) {
        stop(paste0("Considering the name_column_groups '", name_column_groups, "', and that ", as.character(as.character(differential_analysis)), " is passed to differential_analysis and that ", as.character(use_logFC), " is passed to use_logFC,\n",
                    "no suitable Fold Change columns are found in the proteinINFO table"))
      }
    }
  } else {
    if (is.null(Pvalue_column)) {
      stop("Both name_column_groups and Pvalue_column are NULL, so no valid column for Pvalues can be considered from the proteinINFO table")
    } else {
      if (length(Pvalue_column) < 1) {stop("Since name_column_groups is NULL, Pvalue_column must contain at least a valid name")}
      if (any(is.na(Pvalue_column))) {stop("Pvalue_column must not contain NA")}
      if (!is.character(Pvalue_column)) {stop("Pvalue_column must be a character vector with names of columns of proteinINFO")}
      if (!all(Pvalue_column %in% colnames(GCPlist$proteinINFO))) {stop(paste0("Pvalue_column must contain names of columns of proteinINFO. The following names are not present:\n",
                                                                               paste0(Pvalue_column[which(!Pvalue_column%in%colnames(GCPlist$proteinINFO))], collapse = "\n")))}
    }

    if (differential_analysis) {
      if (is.null(FC_column)) {
        stop("You specified differential_analysis as TRUE, but both name_column_groups and FC_column are NULL, so no valid column for Fold change analysis can be considered from the proteinINFO table")
      } else {
        if (length(FC_column)<1) {stop("Since you specified differential_analysis as TRUE, and since name_column_groups is NULL, FC_column must contain at least a valid name")}
        if (any(is.na(FC_column))) {stop("FC_column must not contina missing values")}
        if (!is.character(FC_column)) {stop("FC_column must be a character vector with names of columns of proteinINFO")}
        if (!all(FC_column %in% colnames(GCPlist$proteinINFO))) {stop(paste0("FC_column must contain names of columns of proteinINFO. The following names are not present:\n",
                                                                             paste0(FC_column[which(!FC_column%in%colnames(GCPlist$proteinINFO))], collapse = "\n")))}
      }
    }
  }

  double_Pvalue_columns <- character()
  for (pvc in Pvalue_column) {
    if (length(which(colnames(GCPlist$proteinINFO)==pvc))!=1) {double_Pvalue_columns <- c(double_Pvalue_columns, pvc)}
  }
  if (length(double_Pvalue_columns)>0) {
    stop(paste0("the name following names passed to Pvalue_column are present for more than one column of proteinINFO!\n",
                paste0(double_Pvalue_columns, collapse = "\n")))
  }

  if (differential_analysis) {
    double_FC_column <- character()
    for (fcc in FC_column) {
      if (length(which(colnames(GCPlist$proteinINFO)==fcc))!=1) {double_FC_column <- c(double_FC_column, fcc)}
    }
    if (length(double_FC_column)>0) {
      stop(paste0("the name following names passed to FC_column are present for more than one column of proteinINFO!\n",
                  paste0(double_FC_column, collapse = "\n")))
    }
  }


  if (!all(unname(map_lgl(GCPlist$proteinINFO[, Pvalue_column], is.numeric)))) {stop(paste0("The columns of proteinINFO considered as Pvalue_column must be numeric! The following are not:\n",
                                                                                            paste0(colnames(GCPlist$proteinINFO[, Pvalue_column])[which(!unname(map_lgl(GCPlist$proteinINFO[, Pvalue_column], is.numeric)))], collapse = "\n")))}
  Pvalue_columns_with_weird_num <- character()
  for (this_pvc in Pvalue_column) {
    these_num_no_NA <- pull(GCPlist$proteinINFO, this_pvc)[which(!is.na(pull(GCPlist$proteinINFO, this_pvc)))]
    if (any(these_num_no_NA > 1)  | any(these_num_no_NA < 0)) {
      Pvalue_columns_with_weird_num <- c(Pvalue_columns_with_weird_num, this_pvc)
    }
  }
  if (length(Pvalue_columns_with_weird_num)>0) {
    stop(paste0("The following Pvalue columns contatin numbers that are grather than 1 or less than 0\n",
                paste0(Pvalue_columns_with_weird_num, collapse = "\n")))
  }


  if (length(Pvalue_cutoff)!=1) {stop("Pvalue_cutoff must be a number")}
  if (is.na(Pvalue_cutoff)) {stop("Pvalue_cutoff must be a number")}
  if (!is.numeric(Pvalue_cutoff)) {stop("Pvalue_cutoff must be a number")}
  if (Pvalue_cutoff > 1 | Pvalue_cutoff < 0) {stop("Pvalue_cutoff must be a number between 0 and 1")}

  if (differential_analysis) {
    if (length(FC_cutoff)!=1) {stop("FC_cutoff must be a number")}
    if (is.na(FC_cutoff)) {stop("FC_cutoff must be a number")}
    if (!is.numeric(FC_cutoff)) {stop("FC_cutoff must be a number")}
    if (!is.finite(FC_cutoff)) {stop("FC_cutoff must be finite")}

    if (!all(unname(map_lgl(GCPlist$proteinINFO[, FC_column], is.numeric)))) {stop(paste0("The columns of proteinINFO considered as FC_column must be numeric! The following are not:\n",
                                                                                          paste0(colnames(GCPlist$proteinINFO[, FC_column])[which(!unname(map_lgl(GCPlist$proteinINFO[, FC_column], is.numeric)))], collapse = "\n")))}

    FC_columns_with_infinite_num <- character()
    for (this_fcccc in FC_column) {
      these_num_noNA <- pull(GCPlist$proteinINFO, this_fcccc)[which(!is.na(pull(GCPlist$proteinINFO, this_fcccc)))]
      if (any(!is.finite(these_num_noNA))) {
        FC_columns_with_infinite_num <- c(FC_columns_with_infinite_num, this_fcccc)
      }
    }
    if (length(FC_columns_with_infinite_num)>0) {
      stop(paste0("The following FC columns contatin infinite numbers:\n",
                  paste0(FC_columns_with_infinite_num, collapse = "\n")))
    }
  }

  if (length(background_from_data)!=1) {stop("background_from_data must be TRUE or FALSE")}
  if (is.na(background_from_data)) {stop("background_from_data must be TRUE or FALSE")}
  if (!is.logical(background_from_data)) {stop("background_from_data must be TRUE or FALSE")}

  if (!inherits(OrgDb, "OrgDb")) {
    if (!is.character(OrgDb)) {stop('OrgDb must be an OrgDb object or one of: "human", "mouse", "rat"')}
    if (!identical(tolower(OrgDb), c("human", "mouse", "rat"))) {
      if (length(OrgDb) != 1) {stop('OrgDb must be an OrgDb object or one of: "human", "mouse", "rat"')}
      if (is.na(OrgDb)) {stop('OrgDb must be an OrgDb object or one of: "human", "mouse", "rat"')}
    }
    OrgDb <- tolower(OrgDb)
    OrgDb <- match.arg(OrgDb, c("human", "mouse", "rat"))

    if (OrgDb == "human") {
      needed_package <- "org.Hs.eg.db"
    } else if (OrgDb == "mouse") {
      needed_package <- "org.Mm.eg.db"
    } else if (OrgDb == "rat") {
      needed_package <- "org.Rn.eg.db"
    }

    if (!requireNamespace(needed_package, quietly = TRUE)) {
      stop(paste0('Package "', needed_package, '" is required for OrgDb = "', OrgDb, '". Please, first install it with\n BiocManager::install("', needed_package, '")'))
    }

    OrgDb <- getExportedValue(needed_package, needed_package)
  }

  if (length(keyType)!=1) {stop("keyType must be a character element")}
  if (is.na(keyType)) {stop("keyType must be a character element")}
  if (!is.character(keyType)) {stop("keyType must be a character element")}
  if (!keyType %in% AnnotationDbi::keytypes(OrgDb)) {stop("keyType must be one of the key types supported by OrgDb")}

  valid_symbols <- AnnotationDbi::keys(OrgDb, keytype = keyType)

  if (!is.character(ont)) {stop('ont must be one of: "ALL", "BP", "MF", "CC"')}
  if (!identical(toupper(ont), c("ALL", "BP", "MF", "CC"))) {
    if (length(ont) != 1) {stop('ont must be one of: "ALL", "BP", "MF", "CC"')}
    if (is.na(ont)) {stop('OrgDb must be one of: "ALL", "BP", "MF", "CC"')}
  }
  ont <- toupper(ont)
  ont <- match.arg(ont, c("ALL", "BP", "MF", "CC"))

  if (length(pvalueCutoff)!=1) {stop("pvalueCutoff must be a number")}
  if (is.na(pvalueCutoff)) {stop("pvalueCutoff must be a number")}
  if (!is.numeric(pvalueCutoff)) {stop("pvalueCutoff must be a number")}
  if (pvalueCutoff > 1 | pvalueCutoff < 0) {stop("pvalueCutoff must be a number between 0 and 1")}

  if (!is.character(pAdjustMethod)) {stop('pAdjustMethod must be one of: "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"')}
  if (!identical(pAdjustMethod, c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))) {
    if (length(pAdjustMethod) != 1) {stop('pAdjustMethod must be one of: "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"')}
    if (is.na(pAdjustMethod)) {stop('OrgDb must be one of: "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"')}
  }
  pAdjustMethod <- match.arg(pAdjustMethod, c("none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr"))

  if (length(qvalueCutoff)!=1) {stop("qvalueCutoff must be a number")}
  if (is.na(qvalueCutoff)) {stop("qvalueCutoff must be a number")}
  if (!is.numeric(qvalueCutoff)) {stop("qvalueCutoff must be a number")}
  if (qvalueCutoff > 1 | qvalueCutoff < 0) {stop("qvalueCutoff must be a number between 0 and 1")}



  if (differential_analysis) {

    Pvalue_column_vss <- stringr::str_count(Pvalue_column, "_vs_")
    names(Pvalue_column_vss) <- Pvalue_column
    if (any(Pvalue_column_vss!=1)) {stop(paste0("the column names passed considered for Pvalue_column must have exactly one '_vs_' pattern in their names. The following are not suitable:\n",
                                                paste0(names(Pvalue_column_vss)[which(unname(Pvalue_column_vss!=1))], collapse = "\n")))}

    Pvalue_column_tomatch <- Pvalue_column %>%
      str_remove_all("ttest_") %>%
      str_remove_all("ANOVA_") %>%
      str_remove_all("_PvaluesFDR") %>%
      str_remove_all("_PvalueFDR") %>%
      str_remove_all("_Pvalues") %>%
      str_remove_all("_Pvalue")

    if (any(duplicated(Pvalue_column_tomatch))) {
      Pvalue_column_tomatch_duplicated <- unique(Pvalue_column_tomatch[which(duplicated(Pvalue_column_tomatch))])
      stop(paste0("The following names considered for Pvalue_column produces duplicated:\n",
                  paste0(Pvalue_column[which(Pvalue_column_tomatch%in%Pvalue_column_tomatch_duplicated)])))
    }

    FC_column_vss <- stringr::str_count(FC_column, "_vs_")
    names(FC_column_vss) <- FC_column
    if (any(FC_column_vss!=1)) {stop(paste0("the column names passed considered for FC_column must have exactly one '_vs_' pattern in their names. The following are not suitable:\n",
                                            paste0(names(FC_column_vss)[which(unname(FC_column_vss!=1))], collapse = "\n")))}

    FC_column_tomatch <- FC_column %>%
      str_remove_all("_FC") %>%
      str_remove_all("_logFC")


    if (any(duplicated(FC_column_tomatch))) {
      FC_column_tomatch_duplicated <- unique(FC_column_tomatch[which(duplicated(FC_column_tomatch))])
      stop(paste0("The following names considered for FC_column produces duplicated:\n",
                  paste0(FC_column[which(FC_column_tomatch%in%FC_column_tomatch_duplicated)])))
    }

    matching_tibble <- tibble(Pvalue_column_index = seq(length(Pvalue_column)),
                              matched_comparison = Pvalue_column_tomatch,
                              FC_column_index = rep(NA_integer_, length(Pvalue_column)))

    FC_col_w_no_match <- integer()
    FC_col_w_multiple_matches <- integer()
    for (fcci in seq(length(FC_column))) {

      this_match_with_Pvalue_column_index <- which(matching_tibble$matched_comparison == FC_column_tomatch[fcci])

      if (length(this_match_with_Pvalue_column_index)==0) {
        FC_col_w_no_match <- c(FC_col_w_no_match, fcci)
      } else if (length(this_match_with_Pvalue_column_index)>1) {
        FC_col_w_multiple_matches <- c(FC_col_w_multiple_matches, fcci)
      } else {
        matching_tibble$FC_column_index[this_match_with_Pvalue_column_index] <- fcci
      }
    }

    if (length(FC_col_w_no_match>0)) {stop(paste0("The following FC_column have no match with Pvalue_column:\n",
                                                  paste0(FC_column[FC_col_w_no_match], collapse = "\n")))}
    if (length(FC_col_w_multiple_matches>0)) {stop(paste0("The following FC_column have more than one match with Pvalue_column:\n",
                                                          paste0(FC_column[FC_col_w_multiple_matches], collapse = "\n")))}

    if (anyNA(matching_tibble$FC_column_index)) {
      stop(paste0("The following comparisons have no matching FC column:\n",
                  paste(matching_tibble$matched_comparison[is.na(matching_tibble$FC_column_index)], collapse = "\n")))
    }

    cat("\n\nThe Pvalue_column and the FC_column considered are the following, and matched in this way:\n ")
    matching_tibble_names <- tibble(Space = rep(" ", nrow(matching_tibble)),
                                    Pvalue_column = Pvalue_column[matching_tibble$Pvalue_column_index],
                                    Arrow = rep(" <-->", nrow(matching_tibble)),
                                    FC_column = FC_column[matching_tibble$FC_column_index])
    cat(apply(matching_tibble_names, 1, paste, collapse = "\t"), sep = "\n")

    GO_list_results <- vector(mode = "list", length = 2*nrow(matching_tibble))
    names(GO_list_results)[seq(from = 1, to = length(GO_list_results)-1, by = 2)] <- paste0("higher_in_", matching_tibble$matched_comparison)
    names(GO_list_results)[seq(from = 2, to = length(GO_list_results), by = 2)] <- paste0("lower_in_", matching_tibble$matched_comparison)

  } else {
    matching_tibble <- tibble(Pvalue_column_index = seq(length(Pvalue_column)),
                              matched_comparison = Pvalue_column)
    GO_list_results <- vector(mode = "list", length = nrow(matching_tibble))
    names(GO_list_results) <- matching_tibble$matched_comparison

    if (length(Pvalue_column)==1) {
      cat(paste0("\n\nThe Pvalue_column considered is:\n ", Pvalue_column))
    } else {
      cat(paste0("\n\nThe Pvalue_column considered are:\n",
                 paste0(Pvalue_column, collapse = "\n ")))
    }
  }







  GO_lists_with_all_genes <- vector(mode = "list", length = nrow(matching_tibble))
  names(GO_lists_with_all_genes) <- matching_tibble$matched_comparison






  cat(paste0("\n\nP-values lower than ", Pvalue_cutoff, " are considered significant.\n"))
  if (differential_analysis) {cat(paste0("\nFold Changes > ", FC_cutoff, " are considered for the 'up' GO.\nFold Changes < ", FC_cutoff, " are considered for the 'down' GO.\n"))}
  cat(paste0("\nThe OrgDb ", OrgDb$packageName))
  if (OrgDb$packageName == "org.Hs.eg.db") {cat(" (human)")}
  if (OrgDb$packageName == "org.Mm.eg.db") {cat(" (mouse)")}
  if (OrgDb$packageName == "org.Rn.eg.db") {cat(" (rat)")}
  cat(paste0(" is used; keyType is ", keyType, "; ont is ", ont, ".\n"))
  if (!background_from_data) {cat("Since you set background_from_data as FALSE, all genes listed in the database are used as background!\n")}
  cat(paste0("The pvalueCutoff is ", pvalueCutoff, "; the pAdjustMethod is ", pAdjustMethod, "; the qvalueCutoff is ", qvalueCutoff, ".\n"))



  summary_tibble <- tibble(origin = c("protein_groups_notNA", "unique_gene_names", "valid_gene_names", "signif_protein_groups_notNA", "signif_unique_gene_names", "signif_valid_gene_names"))
  if (differential_analysis) {
    summary_tibble <- add_row(summary_tibble,
                              origin = c("signif_up_protein_groups_notNA", "signif_up_unique_gene_names", "signif_up_valid_gene_names", "signif_down_protein_groups_notNA", "signif_down_unique_gene_names", "signif_down_valid_gene_names"))
  }

  for (tpvci in matching_tibble$Pvalue_column_index) {

    cat(paste0("\n\n____\nGene Ontologies ", tpvci, " of ", length(matching_tibble$Pvalue_column_index), ":  ", matching_tibble$matched_comparison[tpvci], "\n"))

    this_list_with_all_genes <- list(protein_groups_notNA = character(),
                                     unique_gene_names = character(),
                                     valid_gene_names = character(),
                                     signif_protein_groups_notNA = character(),
                                     signif_unique_gene_names = character(),
                                     signif_valid_gene_names = character())

    proteinINFO_rows_to_keep <- !is.na(pull(GCPlist$proteinINFO, Pvalue_column[tpvci]))

    if (differential_analysis) {
      proteinINFO_rows_to_keep <- proteinINFO_rows_to_keep & !is.na(pull(GCPlist$proteinINFO, FC_column[matching_tibble$FC_column_index[tpvci]]))
    }

    proteinINFO_noNA <- GCPlist$proteinINFO[which(proteinINFO_rows_to_keep), ]

    this_list_with_all_genes$protein_groups_notNA <- proteinINFO_noNA$`Gene names`
    names(this_list_with_all_genes$protein_groups_notNA) <- proteinINFO_noNA$protid
    this_list_with_all_genes$unique_gene_names <- this_list_with_all_genes$protein_groups_notNA[which(!duplicated(this_list_with_all_genes$protein_groups_notNA))]
    this_list_with_all_genes$valid_gene_names <- this_list_with_all_genes$unique_gene_names[which(this_list_with_all_genes$unique_gene_names %in% valid_symbols)]

    if (background_from_data) {
      if (length(this_list_with_all_genes$valid_gene_names) > 0) {
        cat(paste0("\n - protein groups not NA: ", length(this_list_with_all_genes$protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$valid_gene_names), "  <- considered as universe in this GO.\n"))
      } else {
        cat(paste0("\n - protein groups not NA: ", length(this_list_with_all_genes$protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$valid_gene_names), "  <- so this GO is not performed!!\n"))
      }
    }

    proteinINFO_signif <- proteinINFO_noNA[which(pull(proteinINFO_noNA, Pvalue_column[tpvci]) < Pvalue_cutoff), ]

    this_list_with_all_genes$signif_protein_groups_notNA <- proteinINFO_signif$`Gene names`
    names(this_list_with_all_genes$signif_protein_groups_notNA) <- proteinINFO_signif$protid
    this_list_with_all_genes$signif_unique_gene_names <- this_list_with_all_genes$signif_protein_groups_notNA[which(!duplicated(this_list_with_all_genes$signif_protein_groups_notNA))]
    this_list_with_all_genes$signif_valid_gene_names <- this_list_with_all_genes$signif_unique_gene_names[which(this_list_with_all_genes$signif_unique_gene_names %in% valid_symbols)]




    if (!differential_analysis) {


      if (background_from_data) {


        if (length(this_list_with_all_genes$valid_gene_names) > 0) {

          if (length(this_list_with_all_genes$signif_valid_gene_names) > 0) {

            cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "  <- considered in this GO.\n"))

            enrich_results <- suppressMessages(clusterProfiler::enrichGO(gene = this_list_with_all_genes$signif_valid_gene_names,
                                                                         OrgDb = OrgDb,
                                                                         keyType = keyType,
                                                                         ont = ont,
                                                                         pvalueCutoff = pvalueCutoff,
                                                                         pAdjustMethod = pAdjustMethod,
                                                                         universe = this_list_with_all_genes$valid_gene_names,
                                                                         qvalueCutoff = qvalueCutoff,
                                                                         ...))
          } else {

            cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "  <- so this GO is not performed!!\n"))

            enrich_results <- NULL
          }
        } else {

          enrich_results <- NULL

        }

      } else {



        if (length(this_list_with_all_genes$signif_valid_gene_names) > 0) {

          cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "  <- considered in this GO.\n"))

          enrich_results <- suppressMessages(clusterProfiler::enrichGO(gene = this_list_with_all_genes$signif_valid_gene_names,
                                                                       OrgDb = OrgDb,
                                                                       keyType = keyType,
                                                                       ont = ont,
                                                                       pvalueCutoff = pvalueCutoff,
                                                                       pAdjustMethod = pAdjustMethod,
                                                                       qvalueCutoff = qvalueCutoff,
                                                                       ...))
        } else {

          cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "  <- so this GO is not performed!!\n"))

          enrich_results <- NULL
        }
      }

      if (!is.null(enrich_results)) {
        GO_list_results[[tpvci]] <- enrich_results@result
      }


      summary_tibble <- bind_cols(summary_tibble, as_tibble(matrix(c(length(this_list_with_all_genes$protein_groups_notNA),
                                                                     length(this_list_with_all_genes$unique_gene_names),
                                                                     length(this_list_with_all_genes$valid_gene_names),
                                                                     length(this_list_with_all_genes$signif_protein_groups_notNA),
                                                                     length(this_list_with_all_genes$signif_unique_gene_names),
                                                                     length(this_list_with_all_genes$signif_valid_gene_names)),
                                                                   nrow = 6,
                                                                   ncol = 1,
                                                                   dimnames = list(NULL, matching_tibble$matched_comparison[tpvci]))))



    } else {


      proteinINFO_signif_up <- proteinINFO_signif[which(pull(proteinINFO_signif, FC_column[matching_tibble$FC_column_index[tpvci]]) > FC_cutoff), ]

      this_list_with_all_genes$signif_up_protein_groups_notNA <- proteinINFO_signif_up$`Gene names`
      names(this_list_with_all_genes$signif_up_protein_groups_notNA) <- proteinINFO_signif_up$protid
      this_list_with_all_genes$signif_up_unique_gene_names <- this_list_with_all_genes$signif_up_protein_groups_notNA[which(!duplicated(this_list_with_all_genes$signif_up_protein_groups_notNA))]
      this_list_with_all_genes$signif_up_valid_gene_names <- this_list_with_all_genes$signif_up_unique_gene_names[which(this_list_with_all_genes$signif_up_unique_gene_names %in% valid_symbols)]



      proteinINFO_signif_down <- proteinINFO_signif[which(pull(proteinINFO_signif, FC_column[matching_tibble$FC_column_index[tpvci]]) < FC_cutoff), ]

      this_list_with_all_genes$signif_down_protein_groups_notNA <- proteinINFO_signif_down$`Gene names`
      names(this_list_with_all_genes$signif_down_protein_groups_notNA) <- proteinINFO_signif_down$protid
      this_list_with_all_genes$signif_down_unique_gene_names <- this_list_with_all_genes$signif_down_protein_groups_notNA[which(!duplicated(this_list_with_all_genes$signif_down_protein_groups_notNA))]
      this_list_with_all_genes$signif_down_valid_gene_names <- this_list_with_all_genes$signif_down_unique_gene_names[which(this_list_with_all_genes$signif_down_unique_gene_names %in% valid_symbols)]


      if (background_from_data) {


        if (length(this_list_with_all_genes$valid_gene_names) > 0) {

          cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "\n"))


          if (length(this_list_with_all_genes$signif_up_valid_gene_names) > 0) {

            cat(paste0("\n   - significantly higher protein groups: ", length(this_list_with_all_genes$signif_up_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_up_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_up_valid_gene_names), "  <- considered in this GO 'up'.\n"))

            enrich_results_up <- suppressMessages(clusterProfiler::enrichGO(gene =  this_list_with_all_genes$signif_up_valid_gene_names,
                                                                            OrgDb = OrgDb,
                                                                            keyType = keyType,
                                                                            ont = ont,
                                                                            pvalueCutoff = pvalueCutoff,
                                                                            pAdjustMethod = pAdjustMethod,
                                                                            universe = this_list_with_all_genes$valid_gene_names,
                                                                            qvalueCutoff = qvalueCutoff,
                                                                            ...))
          } else {

            cat(paste0("\n   - significantly higher protein groups: ", length(this_list_with_all_genes$signif_up_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_up_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_up_valid_gene_names), "  <- so this GO 'up' is not performed!!\n"))

            enrich_results_up <- NULL
          }



          if (length(this_list_with_all_genes$signif_down_valid_gene_names) > 0) {

            cat(paste0("\n   - significantly lower protein groups: ", length(this_list_with_all_genes$signif_down_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_down_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_down_valid_gene_names), "  <- considered in this GO 'down'.\n"))

            enrich_results_down <- suppressMessages(clusterProfiler::enrichGO(gene = this_list_with_all_genes$signif_down_valid_gene_names,
                                                                              OrgDb = OrgDb,
                                                                              keyType = keyType,
                                                                              ont = ont,
                                                                              pvalueCutoff = pvalueCutoff,
                                                                              pAdjustMethod = pAdjustMethod,
                                                                              universe = this_list_with_all_genes$valid_gene_names,
                                                                              qvalueCutoff = qvalueCutoff,
                                                                              ...))

          } else {

            cat(paste0("\n   - significantly lower protein groups: ", length(this_list_with_all_genes$signif_down_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_down_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_down_valid_gene_names), "  <- so this GO 'down' is not performed!!\n"))

            enrich_results_down <- NULL
          }
        } else {

          enrich_results_up <- NULL

          enrich_results_down <- NULL
        }

      } else {

        cat(paste0("\n - significant protein groups: ", length(this_list_with_all_genes$signif_protein_groups_notNA), ",\n   of which unique gene names: ", length(this_list_with_all_genes$signif_unique_gene_names), ",\n   of which valid gene names: ", length(this_list_with_all_genes$signif_valid_gene_names), "\n"))



        if (length(this_list_with_all_genes$signif_up_valid_gene_names) > 0) {

          cat(paste0("\n   - significantly higher protein groups: ", length(this_list_with_all_genes$signif_up_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_up_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_up_valid_gene_names), "  <- considered in this GO 'up'.\n"))

          enrich_results_up <- suppressMessages(clusterProfiler::enrichGO(gene =  this_list_with_all_genes$signif_up_valid_gene_names,
                                                                          OrgDb = OrgDb,
                                                                          keyType = keyType,
                                                                          ont = ont,
                                                                          pvalueCutoff = pvalueCutoff,
                                                                          pAdjustMethod = pAdjustMethod,
                                                                          qvalueCutoff = qvalueCutoff,
                                                                          ...))
        } else {

          cat(paste0("\n   - significantly higher protein groups: ", length(this_list_with_all_genes$signif_up_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_up_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_up_valid_gene_names), "  <- so this GO 'up' is not performed!!\n"))

          enrich_results_up <- NULL
        }



        if (length(this_list_with_all_genes$signif_down_valid_gene_names) > 0) {

          cat(paste0("\n   - significantly lower protein groups: ", length(this_list_with_all_genes$signif_down_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_down_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_down_valid_gene_names), "  <- considered in this GO 'down'.\n"))

          enrich_results_down <- suppressMessages(clusterProfiler::enrichGO(gene = this_list_with_all_genes$signif_down_valid_gene_names,
                                                                            OrgDb = OrgDb,
                                                                            keyType = keyType,
                                                                            ont = ont,
                                                                            pvalueCutoff = pvalueCutoff,
                                                                            pAdjustMethod = pAdjustMethod,
                                                                            qvalueCutoff = qvalueCutoff,
                                                                            ...))

        } else {

          cat(paste0("\n   - significantly lower protein groups: ", length(this_list_with_all_genes$signif_down_protein_groups_notNA), ",\n     of which unique gene names: ", length(this_list_with_all_genes$signif_down_unique_gene_names), ",\n     of which valid gene names: ", length(this_list_with_all_genes$signif_down_valid_gene_names), "  <- so this GO 'down' is not performed!!\n"))

          enrich_results_down <- NULL
        }

      }


      if (!is.null(enrich_results_up)) {
        GO_list_results[[paste0("higher_in_", matching_tibble$matched_comparison[tpvci])]] <- enrich_results_up@result
      }

      if (!is.null(enrich_results_down)) {
        GO_list_results[[paste0("lower_in_", matching_tibble$matched_comparison[tpvci])]] <- enrich_results_down@result
      }


      summary_tibble <- bind_cols(summary_tibble, as_tibble(matrix(c(length(this_list_with_all_genes$protein_groups_notNA),
                                                                     length(this_list_with_all_genes$unique_gene_names),
                                                                     length(this_list_with_all_genes$valid_gene_names),
                                                                     length(this_list_with_all_genes$signif_protein_groups_notNA),
                                                                     length(this_list_with_all_genes$signif_unique_gene_names),
                                                                     length(this_list_with_all_genes$signif_valid_gene_names),
                                                                     length(this_list_with_all_genes$signif_up_protein_groups_notNA),
                                                                     length(this_list_with_all_genes$signif_up_unique_gene_names),
                                                                     length(this_list_with_all_genes$signif_up_valid_gene_names),
                                                                     length(this_list_with_all_genes$signif_down_protein_groups_notNA),
                                                                     length(this_list_with_all_genes$signif_down_unique_gene_names),
                                                                     length(this_list_with_all_genes$signif_down_valid_gene_names)),
                                                                   nrow = 12,
                                                                   ncol = 1,
                                                                   dimnames = list(NULL, matching_tibble$matched_comparison[tpvci]))))




    }


    GO_lists_with_all_genes[[tpvci]] <- this_list_with_all_genes


  }


  cat("\n____\n DONE!\n\n")

  GO_final_list <- list(considered_gene_names = GO_lists_with_all_genes,
                        SUMMARY_gene_names = summary_tibble,
                        GO_results = GO_list_results)


  return(GO_final_list)
}
