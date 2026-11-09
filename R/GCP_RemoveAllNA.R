#' Remove samples with all missing values intensities
#'
#' Starting from a GCPlist, it removes proteins which intensity is missing in all the samples.
#'
#' @param GCPlist a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or ImportOutputProtDiscov().
#'
#' @return the GCPlist in which each table has potentially a reduced number of rows.
#'
#' @examples
#' \dontrun{
#'
#' GCPlist03 <- GCP_RemoveAllNA(GCPlist02)
#'
#' }
#'
#'
#' @export
GCP_RemoveAllNA <- function(GCPlist) {

  checkGCPlist(GCPlist)


  protid_total <- GCPlist$proteinINFO$protid
  protid_NA <- character()


  for (i in 1:length(protid_total)) {
    if (all(is.na(GCPlist$intensities[i, -1]))) {
      protid_NA <- c(protid_NA, GCPlist$intensities$protid[i])
    }
  }

  cat(paste0("\n___\nThere are ", length(protid_NA), " rows with all NAs in the intensity table, which have been removed.\n\n"))

  GCPoutput <- GCPlist
  GCPoutput$proteinINFO <- GCPoutput$proteinINFO[which(!GCPoutput$proteinINFO$protid%in%protid_NA),]
  GCPoutput$intensities <- GCPoutput$intensities[which(!GCPoutput$intensities$protid%in%protid_NA),]

  cat(paste0('Thus, ', length(GCPoutput$proteinINFO$protid), ' out of ', length(GCPlist$proteinINFO$protid) , ' proteins were kept.\n___\n'))

  return(GCPoutput)
}
