#' Checking if the GCPlist is suitable
#'
#' Throw error if GCPlist has something wrong.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#'
#' @return if everything is ok, it does nothing. It only throws an error in case of something wrong.
#'
#' @export
checkGCPlist <- function(GCPlist) {

  if (!is.list(GCPlist)) {stop("GCPlist must be a list with 3 data frames")}
  if (length(GCPlist) != 3) {stop("GCPlist must be a list with 3 data frames")}
  if (!is.data.frame(GCPlist[[1]])) {stop("GCPlist must be a list with 3 data frames")}
  if (!is.data.frame(GCPlist[[2]])) {stop("GCPlist must be a list with 3 data frames")}
  if (!is.data.frame(GCPlist[[3]])) {stop("GCPlist must be a list with 3 data frames")}
  if (any(duplicated(names(GCPlist)))) {stop('The names of the data frames of GCPlist must be "intensities", "proteinINFO", "sampleINFO"')}
  if (!(all(names(GCPlist) %in% c("intensities", "proteinINFO", "sampleINFO")))) {stop('The names of the data frames of GCPlist must be "proteinINFO", "intensities","sampleINFO"')}
  if (length(which(colnames(GCPlist[["proteinINFO"]]) == "protid")) != 1) {stop('The "proteinINFO" dataframe of GCPlist must have exactly one "protid" column')}
  if (length(which(colnames(GCPlist[["intensities"]]) == "protid")) != 1) {stop('The "intensities" dataframe of GCPlist must have exactly one "protid" column')}
  if (any(duplicated(GCPlist[["proteinINFO"]]$protid))) {stop("There must not be duplicated in the protid of the proteinINFO table")}
  if (any(duplicated(GCPlist[["intensities"]]$protid))) {stop("There must not be duplicated in the protid of the intensities table")}

  if (!(all(GCPlist[["intensities"]]$protid %in% GCPlist[["proteinINFO"]]$protid))) {stop('The protid in intensities must be present in the protid of the proteinINFO table')}

  if (any(duplicated(colnames(GCPlist[["intensities"]])[-which(colnames(GCPlist[["intensities"]])=="protid")]))) {stop("There must not be dupicated in the column names of the intensities data frame in the GCPlist")}
  if (any(duplicated(pull(GCPlist[["sampleINFO"]], 1)))) {stop("There must not be duplicated in the sample names reported in the first column of the sampleINFO data frame of the GCPlist")}

  if (!all(pull(GCPlist[["sampleINFO"]], 1) %in% colnames(GCPlist[["intensities"]])[-which(colnames(GCPlist[["intensities"]])=="protid")])) {stop("The names contained in the first column of the sampleINFO data frame of the GCPlist must correspond to column names in the intensities data frames")}

  if (!all(map_lgl(GCPlist[["intensities"]][,colnames(GCPlist[["intensities"]])[which(colnames(GCPlist[["intensities"]])!="protid")]], is.numeric))) {stop("except for protid, all columns of intensities must contain numeric data")}
}
