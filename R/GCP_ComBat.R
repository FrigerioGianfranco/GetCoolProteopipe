#' Adjust for batch effects
#'
#' It performs adjustment for batch effects using an empirical Bayes framework applying the function ComBat from the sva package.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#' @param batch character of length 1 OR a numeric/factor vector. The name of the column of the sampleINFO table containing the batch indications OR a numeric or factor vector containing the batch indications, that will also be added as column in the sampleINFO data frame under a new 'batch' column.
#' @param ... Additional arguments passed to sva::ComBat.
#'
#' @return The GCPlist with the desired intensity values adjusted for batch effects.
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' # passing a numeric/factor vector directly to the 'batch' argument:
#'
#' GCPlist14b1 <- GCP_ComBat(GCPlist = GCPlist14b,
#'                           batch = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1))
#'
#' # passing a column name of sampleINFO, containing the numeric/factor vector, to the 'batch' argument:
#'
#' GCPlist14b$sampleINFO[, "batches"] <- factor(c("batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1"), levels = c("batch1", "batch2", "batch3"))
#'
#' GCPlist14b2 <- GCP_ComBat(GCPlist = GCPlist14b,
#'                           batch = "batches")
#'
#' }
#'
#'
#' @importFrom sva ComBat
#'
#'
#' @export
GCP_ComBat <- function(GCPlist, batch = "Batch", ...) {

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


  the_matrix <- as.matrix(GCPlist$intensities[, which(colnames(GCPlist$intensities)!="protid")])

  combat_matrix <- ComBat(dat = the_matrix, batch = batch_indications, ...)

  GCPoutput$intensities <- as_tibble(combat_matrix)
  GCPoutput$intensities <- add_column(GCPoutput$intensities,
                                      protid = GCPlist$intensities$protid,
                                      .before = which(colnames(GCPlist$intensities)=="protid"))


  return(GCPoutput)
}
