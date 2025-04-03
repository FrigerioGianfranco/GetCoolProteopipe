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

  if (!is.list(GCPlist)) {stop("GCPlist must be a list with 4 data frames")}
  if (length(GCPlist) != 4) {stop("GCPlist must be a list with 4 data frames")}
  if (!is.data.frame(GCPlist[[1]])) {stop("GCPlist must be a list with 4 data frames")}
  if (!is.data.frame(GCPlist[[2]])) {stop("GCPlist must be a list with 4 data frames")}
  if (!is.data.frame(GCPlist[[3]])) {stop("GCPlist must be a list with 4 data frames")}
  if (!is.data.frame(GCPlist[[4]])) {stop("GCPlist must be a list with 4 data frames")}
  if (any(duplicated(names(GCPlist)))) {stop('The names of the data frames of GCPlist must be "quant_raw", "quant_LFQ", "proteinINFO", "sampleINFO"')}
  if (!(all(names(GCPlist) %in% c("quant_raw", "quant_LFQ", "proteinINFO", "sampleINFO")))) {stop('The names of the data frames of GCPlist must be "proteinINFO", "quant_raw", "quant_LFQ","sampleINFO"')}
  if (length(which(colnames(GCPlist[["proteinINFO"]]) == "protid")) != 1) {stop('The "proteinINFO" dataframe of GCPlist must have exactly one "protid" column')}
  if (length(which(colnames(GCPlist[["quant_raw"]]) == "protid")) != 1) {stop('The "quant_raw" dataframe of GCPlist must have exactly one "protid" column')}
  if (length(which(colnames(GCPlist[["quant_LFQ"]]) == "protid")) != 1) {stop('The "quant_LFQ" dataframe of GCPlist must have exactly one "protid" column')}
  if (any(duplicated(GCPlist[["proteinINFO"]]$protid))) {stop("There must not be duplicated in the protid of the proteinINFO table")}
  if (any(duplicated(GCPlist[["quant_raw"]]$protid))) {stop("There must not be duplicated in the protid of the quant_raw table")}
  if (any(duplicated(GCPlist[["quant_LFQ"]]$protid))) {stop("There must not be duplicated in the protid of the quant_LFQ table")}

  if (!(all(GCPlist[["quant_raw"]]$protid %in% GCPlist[["proteinINFO"]]$protid))) {stop('The protid in quant_raw must be present in the protid of the proteinINFO table')}
  if (!(all(GCPlist[["quant_LFQ"]]$protid %in% GCPlist[["proteinINFO"]]$protid))) {stop('The protid in quant_LFQ must be present in the protid of the proteinINFO table')}

  if (!(all(colnames(GCPlist[["quant_raw"]])[-which(colnames(GCPlist[["quant_raw"]])=="protid")] %in% colnames(GCPlist[["quant_LFQ"]])[-which(colnames(GCPlist[["quant_LFQ"]])=="protid")]) & all(colnames(GCPlist[["quant_LFQ"]])[-which(colnames(GCPlist[["quant_LFQ"]])=="protid")] %in% colnames(GCPlist[["quant_raw"]])[-which(colnames(GCPlist[["quant_raw"]])=="protid")]) & length(colnames(GCPlist[["quant_LFQ"]])[-which(colnames(GCPlist[["quant_LFQ"]])=="protid")]) == length(colnames(GCPlist[["quant_raw"]])[-which(colnames(GCPlist[["quant_raw"]])=="protid")]))) {
    stop("There are some differences among raw samples names and LFQ samples names...!")
  }

  if (any(duplicated(colnames(GCPlist[["quant_raw"]])[-which(colnames(GCPlist[["quant_raw"]])=="protid")]))) {stop("There must not be dupicated in the column names of the quant_raw data frame in the GCPlist")}
  if (any(duplicated(colnames(GCPlist[["quant_LFQ"]])[-which(colnames(GCPlist[["quant_LFQ"]])=="protid")]))) {stop("There must not be dupicated in the column names of the quant_LFQ data frame in the GCPlist")}
  if (any(duplicated(pull(GCPlist[["sampleINFO"]], 1)))) {stop("There must not be duplicated in the sample names reported in the first column of the sampleINFO data frame of the GCPlist")}

  if (!all(pull(GCPlist[["sampleINFO"]], 1) %in% colnames(GCPlist[["quant_raw"]])[-which(colnames(GCPlist[["quant_raw"]])=="protid")])) {stop("The names contained in the first column of the sampleINFO data frame of the GCPlist must correspond to column names in the quant_raw data frames")}
  if (!all(pull(GCPlist[["sampleINFO"]], 1) %in% colnames(GCPlist[["quant_LFQ"]])[-which(colnames(GCPlist[["quant_LFQ"]])=="protid")])) {stop("The names contained in the first column of the sampleINFO data frame of the GCPlist must correspond to column names in the quant_LFQ data frames")}

  if (!all(map_lgl(GCPlist[["quant_raw"]][,colnames(GCPlist[["quant_raw"]])[which(colnames(GCPlist[["quant_raw"]])!="protid")]], is.numeric))) {stop("except for protid, all columns of quant_raw must contain numeric data")}
  if (!all(map_lgl(GCPlist[["quant_LFQ"]][,colnames(GCPlist[["quant_LFQ"]])[which(colnames(GCPlist[["quant_LFQ"]])!="protid")]], is.numeric))) {stop("except for protid, all columns of quant_LFQ must contain numeric data")}

}
