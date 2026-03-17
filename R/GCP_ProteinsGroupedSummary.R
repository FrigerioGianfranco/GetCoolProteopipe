#' Get a summary out after the GCP_ProteinsGrouped function.
#'
#' It creates a small table with the summary of the results of the GCP_ProteinsGrouped function. PLEASE NOTE that you must run the function GCP_ProteinsGrouped before!
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param digits_perc NULL or numeric integer of length 1. The number of digits to round the percentages.
#' @param add_perc_symbol logical. If TRUE, it will add the ' %' symbol at the percentages.
#'
#' @return A tibble with summary of the presence of proteins in indicated groups.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' SummaryTable_ProteinsGrouped_filteringNA <- GCP_ProteinsGroupedSummary(GCPlist08)
#' export_the_table(SummaryTable_ProteinsGrouped_filteringNA)
#'
#' }
#'
#'
#'
#' @export
GCP_ProteinsGroupedSummary <- function(GCPlist, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), digits_perc = NULL, add_perc_symbol = FALSE) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
    cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

    if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) < 2) {stop("The sample groups must be at least 2")}

  } else {
    cat("\n -- The name_column_groups considered is NULL --\n\n")
    stop("please, specifiy the name of the sample group column in the name_column_groups argoument")
  }

  if(!is.null(digits_perc)) {
    if (length(digits_perc)!=1) {stop("digits_perc must be a numeric of length 1")}
    if (!is.numeric(digits_perc)) {stop("digits_perc must be a numeric of length 1")}
    if (is.na(digits_perc)) {stop("digits_perc must be a numeric of length 1, not a missing value")}
    is_integerish <- function(x) {is.finite(x) & x == floor(x)}
    if (!is_integerish(digits_perc)) {stop("digits_perc must be an integer number")}
  }

  if (!is.logical(add_perc_symbol)) {stop("add_perc_symbol must be either TRUE or FALSE")}
  if (length(add_perc_symbol) != 1) {stop("add_perc_symbol must be either TRUE or FALSE")}
  if (is.na(add_perc_symbol)) {stop("add_perc_symbol must be either TRUE or FALSE")}


  present_cols <- paste0("present_", levels(pull(GCPlist$sampleINFO, name_column_groups)))
  combination_col <- paste0("combination_", name_column_groups)

  present_combination_cols <- c(present_cols, combination_col)


  if (!all(present_combination_cols %in% colnames(GCPlist$proteinINFO))) {
    stop(paste0('\nAre you sure you run the function GCP_ProteinsGrouped before??\nBecasue the following columns are supposed to be present in proteinINFO, but there are not!\n "',
                paste0(present_combination_cols[which(!(present_combination_cols %in% colnames(GCPlist$proteinINFO)))], collapse = '", "'), '"'))
  }

  if (!all(map_lgl(GCPlist$proteinINFO[, present_cols], is.logical))) {
    stop(paste0('\nAre you sure you run the function GCP_ProteinsGrouped before??\nBecasue the following columns in proteinINFO are supposed to be logical, but they are not!\n "',
                paste0(present_cols[which(!(map_lgl(GCPlist$proteinINFO[, present_cols], is.logical)))], collapse = '", "'), '"'))
  }

  if (!is.factor(pull(GCPlist$proteinINFO, combination_col))) {
    stop(paste0('\nAre you sure you run the function GCP_ProteinsGrouped before??\nBecasue the following column in proteinINFO is supposed to be a factor, but it is not!\n "',
                combination_col, '"'))
  }



  the_summary_table <- GCPlist$proteinINFO %>%
    group_by(!!sym(combination_col)) %>%
    summarise(N=n()) %>%
    mutate(perc = N/sum(N)*100)


  if (!is.null(digits_perc)) {
    the_summary_table[, "perc"] <- formatC(the_summary_table$perc, format = "f", digits = digits_perc)
  }

  if (add_perc_symbol) {
    the_summary_table[, "perc"] <- paste0(the_summary_table$perc, "%")
  }


  return(the_summary_table)
}
