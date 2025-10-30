#' Adjust for batch effects
#'
#' It performs adjustment for batch effects using an empirical Bayes framework applying the function ComBat from the sva package.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param batch character of length 1 OR a numeric/factor vector. The name of the column of the sampleINFO table containing the batch indications OR a numeric or factor vector containing the batch indications, that will also be added as column in the sampleINFO data frame under a new 'batch' column.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The adjustment for batch effects will be performed only in the specified data intensities.
#' @param ... Additional arguments passed to sva::ComBat.
#'
#' @return The GCPlist with the desired intensity values adjusted for batch effects.
#'
#' @importFrom sva ComBat
#'
#' @export
GCP_ComBat <- function(GCPlist, batch = "Batch", raw_or_LFQ = getOption("GetCoolProteopipe.raw_or_LFQ"), ...) {

  checkGCPlist(GCPlist)

  if (is.null(batch)) {stop("batch must not be NULL")}
  if (!is.character(batch) & !is.numeric(batch) & !is.factor(batch)) {stop("batch must be a character vector with the name of the column of the sampleINFO table containing the batch indications OR a numeric or factor vector containing the batch indications")}
  if (length(batch)<1) {stop("batch must have at least one valid value")}
  if (any(is.na(batch))) {stop("batch must not have missing values")}

  if (is.character(batch)) {
    if (length(batch) != 1) {stop("If a character vector, batch must be of length 1")}
    if (!batch%in%colnames(GCPlist$sampleINFO)) {stop("If a character, batch must be a column of sampleINFO")}
    batch_indications <- pull(GCPlist$sampleINFO, batch)
    add_a_column <- FALSE
  } else {
    batch_indications <- batch
    add_a_column <- TRUE
  }

  if (!is.numeric(batch_indications) & !is.factor(batch_indications)) {stop("the batch indications should be a numeric or a factor vector!")}
  if (any(is.na(batch_indications))) {stop("the batch indications must not contain missing values!")}
  if (length(batch_indications) != length(pull(GCPlist$sampleINFO, 1))) {stop("the batch indications must be of the same length of the samples!")}
  if (is.numeric(batch_indications)) {
    if (!all(batch_indications%% 1 == 0)) {stop("The batch indications must not contain non-integer numbers")}
  }
  if (is.factor(batch_indications)) {
    batch_indications <- as.numeric(batch_indications)
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

  cat("\nThe adjustment for batch effects is performed assigning batches to samples in this way:\n Batch Sample\n")
  for (i in 1:length(pull(GCPlist$sampleINFO, 1))) {
    cat(paste0("    ", as.character(batch_indications[i]), "  ", pull(GCPlist$sampleINFO, 1)[i], "\n"))
  }


  GCPoutput <- GCPlist

  if (add_a_column) {
    if (!"Batch"%in%colnames(GCPoutput$sampleINFO)) {
      GCPoutput$sampleINFO[, "Batch"] <- batch_indications
      cat("\nAlso, a column called 'Batch' has been added to sampleINFO, containing those batch indications.\n")
    } else {
      if (identical(GCPoutput$sampleINFO$Batch, batch_indications)) {
        cat("\n")
      } else {
        the_column_name_to_consider <- "Batch"
        for (b in 2:100) {
          if (b==100) {
            cat("\n")
            break
          } else if (b<10) {
            the_number <- paste0("0", as.character(b))
          } else {
            the_number <- as.character(b)
          }
          the_column_name_to_consider <- paste0("Batch", the_number)
          if (the_column_name_to_consider%in%colnames(GCPoutput$sampleINFO)) {
            if (identical(pull(GCPoutput$sampleINFO, the_column_name_to_consider), batch_indications)) {
              cat("\n")
              break
            } else {
              next
            }
          } else {
            GCPoutput$sampleINFO[, the_column_name_to_consider] <- batch_indications
            cat(paste0("\nAlso, a column called '", the_column_name_to_consider, "' has been added to sampleINFO, containing those batch indications.\n"))
            break
          }
        }
      }
    }
  }
  cat("\n")

  if (raw_or_LFQ == "raw") {
    the_raw_matrix <- as.matrix(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")])

    combat_raw_matrix <- ComBat(dat = the_raw_matrix, batch = batch_indications, ...)

    GCPoutput$quant_raw <- as_tibble(combat_raw_matrix)
    GCPoutput$quant_raw <- add_column(GCPoutput$quant_raw,
                                      protid = GCPlist$quant_raw$protid,
                                      .before = which(colnames(GCPlist$quant_raw)=="protid"))
  } else if (raw_or_LFQ == "lfq") {
    the_LFQ_matrix <- as.matrix(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")])

    combat_LFQ_matrix <- ComBat(dat = the_LFQ_matrix, batch = batch_indications, ...)

    GCPoutput$quant_LFQ <- as_tibble(combat_LFQ_matrix)
    GCPoutput$quant_LFQ <- add_column(GCPoutput$quant_LFQ,
                                      protid = GCPlist$quant_LFQ$protid,
                                      .before = which(colnames(GCPlist$quant_LFQ)=="protid"))
  }

  return(GCPoutput)
}
