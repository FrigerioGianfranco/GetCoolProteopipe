#' Imputation of missing values
#'
#' It performs the imputation of missing values, using the same method of the functions scImpute and tImpute from the PhosR package.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#' @param quant_rate numeric, between 0 an 1 (included). Quantification rate per group, considered for the scImpute function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param seed numeric. The value that will be used for set.seed.
#' @param method one of: "both", "scImpute", "tImpute". If both is selected, it performs first the scImpute, and then the tImpute.
#'
#' @return a GCPlist list with the missing values imputed.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist09 <- GCP_NAimputation(GCPlist = GCPlist08,
#'                               quant_rate = 0.5)
#'
#' }
#'
#'
#'
#' @importFrom PhosR PhosphoExperiment scImpute tImpute
#'
#' @export
GCP_NAimputation <- function(GCPlist, quant_rate = 0.5, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), seed = 123, method = c("both", "scImpute", "tImpute")) {

  checkGCPlist(GCPlist)


  if (length(quant_rate)!=1) {stop("quant_rate must be a numeric of length 1")}
  if (!is.numeric(quant_rate)) {stop("quant_rate must be a numeric of length 1")}
  if (is.na(quant_rate)) {stop("quant_rate must be a numeric of length 1, not a missing value")}
  if (quant_rate<0 | quant_rate>1) {stop("quant_rate must be between 0 and 1 (included)")}


  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}


    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }
    if (any(table(pull(GCPlist$sampleINFO, name_column_groups)) == 0)) {
      GCPlist$sampleINFO[,name_column_groups] <- droplevels(pull(GCPlist$sampleINFO, name_column_groups))
    }


  } else {
    cat("\n -- The name_column_groups considered is NULL --\n\n")
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

    if (is.null(name_column_groups)) {

      cat("\nNote: name_column_groups was not specidied, so scimpute is not performed!\n")

      do_scimpute <- FALSE
      do_timpute <- TRUE

    } else if (any(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n()))$N == 1)) {
      cat(paste0('\nNote: the following groups have only one sample, so scimpute is not performed!\n ',
                     paste0(as.character(pull(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n())), name_column_groups))[which(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n()))$N == 1)], collapse = ", ")))
      cat("\n")


      do_scimpute <- FALSE
      do_timpute <- TRUE

    } else {
      do_scimpute <- TRUE
      do_timpute <- TRUE
    }

  } else if (method == "scimpute") {

    if (is.null(name_column_groups)) {
      stop("name_column_groups was not specidied, so scimpute is not performed")
    } else if (any(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n()))$N == 1)) {
      stop(paste0('the following groups have only one sample, so scimpute cannot be done!\n ',
                  paste0(as.character(pull(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n())), name_column_groups))[which(as_tibble(GCPlist$sampleINFO%>%group_by(!!sym(name_column_groups))%>%summarise(N = n()))$N == 1)], collapse = ", ")))
    } else {
      do_scimpute <- TRUE
      do_timpute <- FALSE
    }

  } else if (method == "timpute") {
    do_scimpute <- FALSE
    do_timpute <- TRUE
  }


  check_ppe_object <- function(ppe_object, the_matrix_type) {

    GCPtable <- GCPlist$intensities

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


  print_missing_info <- function(ppe_object, the_matrix_type) {

    ppe_matrix <- ppe_object@assays@data@listData[[the_matrix_type]]

    cat(paste0("\n- ", sum(is.na(ppe_matrix)), " / ", length(ppe_matrix), " (", round(sum(is.na(ppe_matrix))/length(ppe_matrix)*100, digits = 1), "%)"))
  }


  GCPoutput <- GCPlist

  set.seed(seed)

  ppe_intensities <- suppressWarnings(PhosphoExperiment(assays = list(Quantification = as.matrix(GCPlist$intensities[,which(colnames(GCPlist$intensities)!="protid")]))))

  check_ppe_object(ppe_intensities, "Quantification")

  if (!is.null(name_column_groups)) {
    ppe_intensities@colData@listData[["condition"]] <- pull(GCPlist$sampleINFO, name_column_groups)
  }


  cat("\n______\nIn the intensities table, the number of NAs is:")
  print_missing_info(ppe_intensities, "Quantification")
  cat(", initially.")

  if (do_scimpute) {

    ppe_intensities <- suppressWarnings(scImpute(ppe_intensities, quant_rate, pull(GCPlist$sampleINFO, name_column_groups)))

    check_ppe_object(ppe_intensities, "imputed")

    print_missing_info(ppe_intensities, "imputed")
    cat(", after the scImpute.")
  }

  if (do_timpute) {

    if (do_scimpute) {
      ppe_intensities <- tImpute(ppe_intensities, assay = "imputed")
    } else {
      ppe_intensities <- tImpute(ppe_intensities)
    }

    check_ppe_object(ppe_intensities, "imputed")

    print_missing_info(ppe_intensities, "imputed")
    cat(", after the tImpute.")
  }

  GCPoutput$intensities <- as_tibble(ppe_intensities@assays@data@listData[["imputed"]]) %>%
    add_column(protid = GCPlist$intensities$protid,
               .before = 1)

  cat("\n______\n")

  return(GCPoutput)
}

