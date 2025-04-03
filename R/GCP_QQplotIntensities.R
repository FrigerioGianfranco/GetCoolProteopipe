#' Plotting an overall Q-Q plot
#'
#' It crates a Q-Q plot considering all the intensities of protein from all samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ".  The plot will be created only with the specified data intensities.
#' @param Title NULL or character of length 1. The title you want to ad on the top of the plot.
#'
#' @return a ggplot object.
#'
#' @importFrom ggpubr ggqqplot
#'
#' @export
GCP_QQplotIntensities <- function(GCPlist, raw_or_LFQ = c("raw", "LFQ"), Title = "QQ plot intensities") {

  checkGCPlist(GCPlist)

  if (!identical(tolower(raw_or_LFQ), c("raw", "lfq"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("raw", "lfq"))

  if (!is.null(Title)) {
    if (length(Title)!=1) {
      stop("Pass a suitable Title!")
    } else {
      if (!is.character(Title) & !is.numeric(Title) & !is.logical(Title)) {
        stop("Pass a suitable Title!")
      }
    }
  }

  if (raw_or_LFQ == "raw") {
    all_intensities <- GCPlist$quant_raw[, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw) != "protid")]]  %>%
      gather(key="sample", value="value")
  } else if (raw_or_LFQ == "lfq") {
    all_intensities <- GCPlist$quant_LFQ[, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ) != "protid")]]  %>%
      gather(key="sample", value="value")
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }


  the_QQplot <- ggpubr::ggqqplot(all_intensities$value) +
    ggtitle(Title) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust = 0.5))

  return(the_QQplot)

}
