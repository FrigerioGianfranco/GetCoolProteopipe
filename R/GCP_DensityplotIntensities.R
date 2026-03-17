#' Plotting an overall density plot
#'
#' It crates an histogram and density curve considering all the intensities of protein from all samples.
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
#' Fig01_IntensityDistribution_before_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist05,
#'                                                                             Title = "Distribution of intensities, before processing")
#' export_figures(Fig01_IntensityDistribution_before_processing)
#'
#' }
#'
#'
#'
#' @export
GCP_DensityplotIntensities <- function(GCPlist, Title = "Distribution of intensities") {

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

  binwidth_calculated <- 2 * IQR(all_intensities$value, na.rm = TRUE) / length(all_intensities$value)^(1/3)

  the_density_plot <- ggplot(all_intensities, aes(x=value)) +
    geom_histogram(aes(y=after_stat(density)), binwidth = binwidth_calculated, color="black", fill="grey", alpha=0.6) +
    geom_density(color="red", linewidth=1) +
    stat_function(fun = dnorm, args = list(mean=mean(all_intensities$value, na.rm=TRUE), sd=sd(all_intensities$value, na.rm=TRUE)),
                  color="darkgreen", linetype="dashed", linewidth=1) +
    ggtitle(Title) +
    xlab("Theoretical") +
    ylab("Density") +
    theme_bw() +
    theme(plot.title = element_text(face="bold", hjust = 0.5))


  return(the_density_plot)

}
