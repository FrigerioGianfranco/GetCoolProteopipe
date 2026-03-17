#' Plotting an overall Q-Q plot
#'
#' It crates a Q-Q plot considering all the intensities of protein from all samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param Title NULL or character of length 1. The title you want to ad on the top of the plot.
#'
#' @return a ggplot object.
#'
#'
#' @examples
#' \dontrun{
#'
#' Fig02_QQ_plot_before_processing <- GCP_QQplotIntensities(GCPlist = GCPlist05,
#'                                                          Title = "QQ plot - before processing")
#' export_figures(Fig02_QQ_plot_before_processing)
#'
#' }
#'
#'
#'
#' @importFrom ggpubr ggqqplot
#'
#' @export
GCP_QQplotIntensities <- function(GCPlist, Title = "QQ plot intensities") {

  checkGCPlist(GCPlist)

  if (!is.null(Title)) {
    if (length(Title)!=1) {
      stop("Pass a suitable Title!")
    } else {
      if (!is.character(Title) & !is.numeric(Title) & !is.logical(Title)) {
        stop("Pass a suitable Title!")
      }
    }
  }

  all_intensities <- GCPlist$intensities[, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities) != "protid")]]  %>%
    gather(key="sample", value="value")

  the_QQplot <- ggpubr::ggqqplot(all_intensities$value) +
    ggtitle(Title) +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust = 0.5))

  return(the_QQplot)

}
