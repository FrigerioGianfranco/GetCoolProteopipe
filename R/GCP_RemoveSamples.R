#' Remove samples from the dataset
#'
#' Remove samples from both the intensities and  SampleINFO data frames. Only one of the arguments among remove_samples, keep_samples, remove_groups, or keep_groups has to be specified.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param remove_samples NULL or a character containing the name of samples to remove.
#' @param keep_samples NULL or a character containing the name of samples to keep.
#' @param remove_groups NULL or a character containing the name of groups to remove.
#' @param keep_groups NULL or a character containing the name of groups to keep.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups. It must be specified if you used remove_groups or keep_groups.
#'
#' @return the GCPlist with samples removed from the sampleINFO, and intensities data frames.
#'
#'
#' @examples
#' \dontrun{
#'
#' # specify what to remove:
#'
#' GCPlist00rm1 <- GCP_RemoveSamples(GCPlist = GCPlist00,
#'                                   remove_samples = c("S2", "S5"),
#'                                   keep_samples = NULL)
#'
#' # or specify what to keep (same output as above):
#'
#' GCPlist00rm2 <- GCP_RemoveSamples(GCPlist = GCPlist00,
#'                                   remove_samples = NULL,
#'                                   keep_samples = c("S1", "S3", "S4", "V1", "V2", "V3", "V4", "V5"))
#'
#' }
#'
#'
#' @export
GCP_RemoveSamples <- function(GCPlist, remove_samples = NULL, keep_samples = NULL, remove_groups = NULL, keep_groups = NULL, name_column_groups = NULL) {

  checkGCPlist(GCPlist)

  GCPoutput <- GCPlist

  if (is.null(remove_samples) & is.null(keep_samples) & is.null(remove_groups) & is.null(keep_groups)) {

    cat("\n")
    cat("______\n")
    cat("Nothing changed!\n")
    cat("______\n")

    return(GCPoutput)

  } else if (sum(is.null(remove_samples), is.null(keep_samples), is.null(remove_groups), is.null(keep_groups)) != 3) {
    stop("Only one of the arguments among remove_samples, keep_samples, remove_groups, or keep_groups has to be specified")
  } else {

    samples_to_keep <- pull(GCPlist$sampleINFO, 1)
    change_order <- FALSE


    if (!is.null(remove_samples)) {
      if (!is.character(remove_samples)) {stop("remove_samples must be a character vector")}
      if (length(remove_samples)<1) {stop("remove_samples must contain valid names")}
      if (any(is.na(remove_samples))) {stop("remove_samples must not contain missing values")}
      if (!all(remove_samples %in% pull(GCPlist$sampleINFO, 1))) {
        stop(paste0("Names passed to remove_samples must be names of samples. The following sample names you passed are not present:\n",
                    paste0(remove_samples[which(!remove_samples %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
      }

      samples_to_keep <- samples_to_keep[which(!samples_to_keep %in% remove_samples)]

    } else if (!is.null(keep_samples)) {
      if (!is.character(keep_samples)) {stop("keep_samples must be a character vector")}
      if (length(keep_samples)<1) {stop("keep_samples must contain valid names")}
      if (any(is.na(keep_samples))) {stop("keep_samples must not contain missing values")}
      if (!all(keep_samples %in% pull(GCPlist$sampleINFO, 1))) {
        stop(paste0("Names passed to keep_samples must be names of samples. The following sample names you passed are not present:\n",
                    paste0(keep_samples[which(!keep_samples %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
      }

      samples_to_keep <- keep_samples
      change_order <- TRUE

    } else if (!is.null(remove_groups)) {
      if (!is.character(remove_groups)) {stop("remove_groups must be a character vector")}
      if (length(remove_groups)<1) {stop("remove_groups must contain valid names")}
      if (any(is.na(remove_groups))) {stop("remove_groups must not contain missing values")}

      if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
      if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
      if (is.na(name_column_groups)) {stop("name_column_groups must a character of length 1, not a missing value")}
      if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

      if (!all(remove_groups %in% pull(GCPlist$sampleINFO, name_column_groups))) {
        stop(paste0("Names passed to remove_groups must be names of groups in the column ", name_column_groups, " of the sampleINFO dataframe. The following group names you passed are not present:\n",
                    paste0(remove_groups[which(!remove_groups %in% pull(GCPlist$sampleINFO, name_column_groups))], collapse = "\n")))
      }

      samples_to_keep <- pull(GCPlist$sampleINFO, 1)[which(!pull(GCPlist$sampleINFO, name_column_groups) %in% remove_groups)]

    } else if (!is.null(keep_groups)) {
      if (!is.character(keep_groups)) {stop("keep_groups must be a character vector")}
      if (length(keep_groups)<1) {stop("keep_groups must contain valid names")}
      if (any(is.na(keep_groups))) {stop("keep_groups must not contain missing values")}

      if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
      if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
      if (is.na(name_column_groups)) {stop("name_column_groups must a character of length 1, not a missing value")}
      if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

      if (!all(keep_groups %in% pull(GCPlist$sampleINFO, name_column_groups))) {
        stop(paste0("Names passed to keep_groups must be names of groups in the column ", name_column_groups, " of the sampleINFO dataframe. The following group names you passed are not present:\n",
                    paste0(keep_groups[which(!keep_groups %in% pull(GCPlist$sampleINFO, name_column_groups))], collapse = "\n")))
      }

      samples_to_keep <- pull(GCPlist$sampleINFO, 1)[which(pull(GCPlist$sampleINFO, name_column_groups) %in% keep_groups)]

    } else {
      stop("One of the arguments among remove_samples, keep_samples, remove_groups, or keep_groups has to be specified")
    }


    GCPoutput$sampleINFO <- GCPoutput$sampleINFO[which(pull(GCPoutput$sampleINFO, 1) %in% samples_to_keep),]
    GCPoutput$intensities <- GCPoutput$intensities[, c("protid", samples_to_keep)]

    if (change_order) {
      GCPoutput$sampleINFO <- GCPoutput$sampleINFO[match(samples_to_keep, pull(GCPoutput$sampleINFO, 1)),]
    }

    if (all(pull(GCPlist$sampleINFO, 1) %in% pull(GCPoutput$sampleINFO, 1))) {

      cat("\n")
      cat("______\n")
      cat("Nothing changed!\n")
      cat("______\n")

    } else {

      cat("\n")
      cat("______\n")
      cat("- The following samples have been kept:\n")
      cat(paste0(pull(GCPoutput$sampleINFO, 1), collapse = "\n"))
      cat("\n\n")
      cat("- The following samples have been removed:\n")
      cat(paste0(pull(GCPlist$sampleINFO, 1)[which(!pull(GCPlist$sampleINFO, 1) %in% pull(GCPoutput$sampleINFO, 1))], collapse = "\n"))
      cat("\n______\n")

    }


    return(GCPoutput)

  }
}


