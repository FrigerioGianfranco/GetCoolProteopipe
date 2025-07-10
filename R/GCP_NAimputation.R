#' Imputation of missing values
#'
#' It performs the imputation of missing values, using the same method of the functions scImpute and tImpute from the PhosR package. MORE DESCRIPTION!
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.#'
#' @param quant_rate numeric, between 0 an 1 (included). Quantification rate per group, considered for the scImpute function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param seed numeric. The value that will be used for set.seed.
#' @param method one of: "both", "scImpute", "tImpute". If both is selected, it performs first the scImpute, and then the tImpute.
#'
#' @return a GCPlist list with the missing values imputed.
#'
#' @importFrom PhosR PhosphoExperiment scImpute tImpute
#'
#' @export
GCP_NAimputation <- function(GCPlist, quant_rate = 0.5, name_column_groups = NULL, seed = 123, method = c("both", "scImpute", "tImpute")) {

  checkGCPlist(GCPlist)


  if (length(quant_rate)!=1) {stop("quant_rate must be a numeric of length 1")}
  if (!is.numeric(quant_rate)) {stop("quant_rate must be a numeric of length 1")}
  if (is.na(quant_rate)) {stop("quant_rate must be a numeric of length 1, not a missing value")}
  if (quant_rate<0 | quant_rate>1) {stop("quant_rate must be between 0 and 1 (included)")}

  allwiththis_created <- FALSE

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

  } else {
    if ("allwiththis" %in% colnames(GCPlist$sampleINFO)) {stop("Please, don't call a column of the sampleINFO data frame 'allwiththis' as I need to create one with this name now")}
    GCPlist$sampleINFO <- mutate(GCPlist$sampleINFO, allwiththis = as.factor("theOnlyGroup"))

    name_column_groups <- "allwiththis"

    allwiththis_created <- TRUE
  }

  if (length(seed)!=1) {stop("seed must be a numeric of length 1")}
  if (!is.numeric(seed)) {stop("seed must be a numeric of length 1")}
  if (is.na(seed)) {stop("seed must be a numeric of length 1, not a missing value")}

  if (!identical(tolower(method), c("both", "scimpute", "timpute"))) {
    if (length(method) != 1) {stop('method must be one of "both", "scImpute", "tImpute"')}
    if (is.na(method)) {stop('method must be one of "both", "scImpute", "tImpute"')}
  }
  method <- tolower(method)
  method <- match.arg(method, c("both", "scimpute", "timpute"))

  if (method == "both") {
    do_scimpute <- TRUE
    do_timpute <- TRUE
  } else if (method == "scimpute") {
    do_scimpute <- TRUE
    do_timpute <- FALSE
  } else if (method == "timpute") {
    do_scimpute <- FALSE
    do_timpute <- TRUE
  }


  check_ppe_object <- function(ppe_object, rl, the_matrix_type) {

    if (rl == "raw") {GCPtable <- GCPlist$quant_raw} else if (rl == "LFQ") {GCPtable <- GCPlist$quant_LFQ}

    if (!is.matrix(ppe_object@assays@data@listData[[the_matrix_type]])) {
      stop("For some reason there was not a matrix inside the created ppe object... that's wired, ask Gianfranco to check!!")
    }
    if (!identical(colnames(ppe_object@assays@data@listData[[the_matrix_type]]), colnames(GCPtable)[which(colnames(GCPtable)!="protid")])) {
      stop("For some reason the column names were changed while creating the ppe object... that's wired, ask Gianfranco to check!!")
    }
    if (nrow(ppe_object@assays@data@listData[[the_matrix_type]]) != nrow(GCPtable)) {
      stop("For some reason the row numbers were changed while creating the ppe object... that's wired, ask Gianfranco to check!!")
    }
  }


  print_missing_info <- function(ppe_object, rl, the_matrix_type) {

    ppe_matrix <- ppe_object@assays@data@listData[[the_matrix_type]]

    cat(paste0("\n- ", sum(is.na(ppe_matrix)), " / ", length(ppe_matrix), " (", round(sum(is.na(ppe_matrix))/length(ppe_matrix)*100, digits = 1), "%)"))
  }


  GCPoutput <- GCPlist

  set.seed(seed)

  # raw

  ppe_raw <- suppressWarnings(PhosphoExperiment(assays = list(Quantification = as.matrix(GCPlist$quant_raw[,which(colnames(GCPlist$quant_raw)!="protid")]))))

  check_ppe_object(ppe_raw, "raw", "Quantification")

  ppe_raw@colData@listData[["condition"]] <- pull(GCPlist$sampleINFO, name_column_groups)

  cat("\n______\nIn the raw table, the number of NAs is:")
  print_missing_info(ppe_raw, "raw", "Quantification")
  cat(", initially.")

  if (do_scimpute) {

    ppe_raw <- suppressWarnings(scImpute(ppe_raw, quant_rate, pull(GCPlist$sampleINFO, name_column_groups)))

    check_ppe_object(ppe_raw, "raw", "imputed")

    print_missing_info(ppe_raw, "raw", "imputed")
    cat(", after the scImpute.")
  }

  if (do_timpute) {

    if (do_scimpute) {
      ppe_raw <- tImpute(ppe_raw, assay = "imputed")
    } else {
      ppe_raw <- tImpute(ppe_raw)
    }

    check_ppe_object(ppe_raw, "raw", "imputed")

    print_missing_info(ppe_raw, "raw", "imputed")
    cat(", after the tImpute.")
  }

  GCPoutput$quant_raw <- as_tibble(ppe_raw@assays@data@listData[["imputed"]]) %>%
    add_column(protid = GCPlist$quant_raw$protid,
               .before = 1)


  # LFQ

  ppe_LFQ <- suppressWarnings(PhosphoExperiment(assays = list(Quantification = as.matrix(GCPlist$quant_LFQ[,which(colnames(GCPlist$quant_LFQ)!="protid")]))))

  check_ppe_object(ppe_LFQ, "LFQ", "Quantification")

  ppe_LFQ@colData@listData[["condition"]] <- pull(GCPlist$sampleINFO, name_column_groups)

  cat("\n\n____\nIn the LFQ table, the number of NAs is:")
  print_missing_info(ppe_LFQ, "LFQ", "Quantification")
  cat(", initially.")


  if (do_scimpute) {

    ppe_LFQ <- suppressWarnings(scImpute(ppe_LFQ, quant_rate, pull(GCPlist$sampleINFO, name_column_groups)))

    check_ppe_object(ppe_LFQ, "LFQ", "imputed")

    print_missing_info(ppe_LFQ, "LFQ", "imputed")
    cat(", after the scImpute.")
  }

  if (do_timpute) {

    if (do_scimpute) {
      ppe_LFQ <- tImpute(ppe_LFQ, assay = "imputed")
    } else {
      ppe_LFQ <- tImpute(ppe_LFQ)
    }


    check_ppe_object(ppe_LFQ, "LFQ", "imputed")

    print_missing_info(ppe_LFQ, "LFQ", "imputed")
    cat(", after the tImpute.")
  }


  GCPoutput$quant_LFQ <- as_tibble(ppe_LFQ@assays@data@listData[["imputed"]]) %>%
    add_column(protid = GCPlist$quant_LFQ$protid,
               .before = 1)

  cat("\n______\n")

  if (allwiththis_created) {
    GCPoutput$sampleINFO <- GCPoutput$sampleINFO[, which(colnames(GCPoutput$sampleINFO) != "allwiththis")]
  }

  return(GCPoutput)

}

