#' Remove samples with all zero intensities
#'
#' Starting from a GCPlist, it removes proteins which intensity is equal to zero in all the samples.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#'
#' @return the GCPlist in which each table has potentially a reduced number of rows.
#'
#' @examples
#' \dontrun{
#'
#' GCPlist03 <- GCP_RemoveAllZero(GCPlist02)
#'
#' }
#'
#'
#' @export
GCP_RemoveAllZero <- function(GCPlist) {

  checkGCPlist(GCPlist)


  protid_total <- GCPlist$proteinINFO$protid
  protid_zero <- character()


  for (i in 1:length(protid_total)) {
    if (all(as.numeric(GCPlist$intensities[i, -1])==0)) {
      protid_zero <- c(protid_zero, GCPlist$intensities$protid[i])
    }
  }

  cat(paste0("\n___\nThere are ", length(protid_zero), " rows with all zero in the intensity table, which have been removed.\n\n"))

  GCPoutput <- GCPlist
  GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_zero),]
  GCPoutput$intensities <- GCPoutput$intensities[which(!GCPoutput$intensities$protid%in%protid_zero),]

  cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept.\n___\n'))

  return(GCPoutput)
}
