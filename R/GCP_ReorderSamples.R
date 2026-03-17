#' Reorder samples
#'
#' Reorder samples in the CGP list.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param sample_names_ordered NULL or a character or factor vector containing the name of samples with the new desired order.
#' @param name_column_groups NULL or a character of length 1. The name of the column of the sampleINFO table containing the sample groups. If sample_names_ordered is NULL, the order will be based on the on such groups.
#'
#' @return the GCPlist with samples reordered.
#'
#' @examples
#' \dontrun{
#'
#' # Reordering samples based on the groups (by default based on the set name_column_groups):
#'
#' GCPlist02rd1 <- GCP_ReorderSamples(GCPlist02)
#'
#' # Reordering samples specifying the samples:
#'
#' GCPlist02rd2 <- GCP_ReorderSamples(GCPlist02, sample_names_ordered = c("S_1", "V_1", "S_2", "V_2", "S_3", "V_3", "S_4", "V_4", "S_5", "V_5"))
#'
#'
#' }
#'
#'
#' @export
GCP_ReorderSamples <- function(GCPlist, sample_names_ordered = NULL, name_column_groups = getOption("GetCoolProteopipe.name_column_groups")) {

  checkGCPlist(GCPlist)

  if (!is.null(sample_names_ordered)) {
    if (length(sample_names_ordered) != length(pull(GCPlist$sampleINFO, 1))) {stop("sample_names_ordered must have the same length of the previous sample names!")}
    if (any(is.na(sample_names_ordered))) {stop("sample_names_ordered must not contian NAs")}
    if (!is.character(sample_names_ordered) & !is.factor(sample_names_ordered)) {stop("sample_names_ordered must be a character or factor vector")}
    if (any(duplicated(sample_names_ordered))) {stop(paste0("sample_names_ordered must not contain duplicated. The followings are duplicated:\n",
                                                            paste0(unique(sample_names_ordered[which(duplicated(sample_names_ordered))]), collapse = "\n")))}
    if (any(!sample_names_ordered%in%pull(GCPlist$sampleINFO, 1))) {stop(paste0("All elements of sample_names_ordered must be present in the first column of the sampleINFO table. The followings are not:\n",
                                                                                paste0(sample_names_ordered[which(!sample_names_ordered%in%pull(GCPlist$sampleINFO, 1))], collapse = "\n")))}

    GCPoutput <- GCPlist

    GCPoutput$sampleINFO <- arrange(GCPlist$sampleINFO, match(.data[[colnames(GCPlist$sampleINFO)[1]]], sample_names_ordered))

    column_sample_names_ordered <- c("protid", sample_names_ordered)

    GCPoutput$intensities <- select(GCPlist$intensities, all_of(column_sample_names_ordered))

  } else {
    if (is.null(name_column_groups)) {stop("both the arguments sample_names_ordered and name_column_groups are NULL: so no sample reordering can be performed!")}
    if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
    cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

    GCPoutput <- GCPlist

    GCPoutput$sampleINFO <- arrange(GCPlist$sampleINFO, !!sym(name_column_groups))

    column_sample_names_ordered <- c("protid", pull(GCPoutput$sampleINFO, 1))

    GCPoutput$intensities <- select(GCPlist$intensities, all_of(column_sample_names_ordered))

  }

  if (all(pull(GCPoutput$sampleINFO, 1) == pull(GCPlist$sampleINFO, 1))) {
    cat("\n______\nNothing changed!\n______\n")
  } else {
    cat("\n______\nThe samples have been reordered in this way:\n ")
    cat(paste0(pull(GCPoutput$sampleINFO, 1), collapse = "\n "))
    cat("\n______\n")
  }

  return(GCPoutput)
}
