#' Reorder samples
#'
#' Reorder samples in the CGP list.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param sample_names_ordered a character or factor vector containing the name of samples with the new desired order.
#'
#' @return the GCPlist with samples reordered.
#'
#' @export
GCP_ReorderSamples <- function(GCPlist, sample_names_ordered) {
  
  checkGCPlist(GCPlist)
  
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
  
  GCPoutput$quant_raw <- select(GCPlist$quant_raw, all_of(column_sample_names_ordered))
  GCPoutput$quant_LFQ <- select(GCPlist$quant_LFQ, all_of(column_sample_names_ordered))
  
  return(GCPoutput)
  
}
