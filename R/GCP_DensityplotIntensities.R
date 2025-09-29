#' Plotting an overall density plot
#'
#' It crates an histogram and density curve considering all the intensities of protein from all samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The plot will be created only with the specified data intensities.
#' @param Title NULL or character of length 1. The title you want to ad on the top of the plot.
#'
#' @return a ggplot object.
#'
#' @export
GCP_DensityplotIntensities <- function(GCPlist, raw_or_LFQ = getOption("GetCoolProteopipe.raw_or_LFQ"), Title = "Distribution of intensities") {

  checkGCPlist(GCPlist)

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
