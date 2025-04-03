#' Creating a Volcano plot
#'
#' It performs a Volcano plot on the results of a t-test and a Fold Change analysis.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function. T-test and Fold Change analyses must have been performed with the functions GCP_ttest and GCP_FoldChange.
#' @param FDR logical. If TRUE, the plot will be build using the FDR p-values; if FALSE, the p-values.
#' @param log_base numeric of length 1. The log-base passed to GCP_FoldChange.
#' @param pcutoff_colored NULL or a numeric of length 1. The p-value cut-off below which showing colored dots. If NULL, dots will not be colored by significance.
#' @param pcutoff_line NULL or a numeric of length 1. The p-value where a horizontal dahsed line will be added to the plot. If NULL, this line will not be showed.
#' @param pcutoff_prot_label NULL or a numeric of length 1. The p-value below which protein names will be added to the plot. If NULL, no names will be reported.
#' @param name_column_proteinlabels NULL or a character of length 1. The name of a column in the proteinINFO table containing the names you would plot.
#' @param name_column_proteingroups NULL or a character of length 1. The name of a column in the proteinINFO table containing groups of protein. If you pass this argument, the dots will be colored based on these groups, and this will replace whether you specify pcutoff_colored.
#' @param col_pal_difference a character vector containing exactly 3 colors. If pcutoff_colored is specified, they will be used for not significant, decreased, and increased proteins, respectively.
#' @param col_pal_groups NULL or a character vector with colors to be used if name_column_proteingroups were specified. If NULL, colors from the pals package will be used (see function GetFeatistics::build_long_vector_of_colors).
#'
#' @return A ggplot object.
#'
#' @export
GCP_Volcano <- function(GCPlist, FDR = TRUE, log_base = 2, pcutoff_colored = 0.05, pcutoff_line = 0.05, pcutoff_prot_label = NULL, name_column_proteinlabels = NULL, name_column_proteingroups = NULL, col_pal_difference = c("grey", "blue", "red"), col_pal_groups = NULL) {

  checkGCPlist(GCPlist)

  if (!"logFC" %in% colnames(GCPlist$proteinINFO)) {stop('"logFC" is not present among the columns of proteinINFO. Did you run GCP_FoldChange?')}
  if (!"FCcomparison" %in% colnames(GCPlist$proteinINFO)) {stop('"FCcomparison" is not present among the columns of proteinINFO. Did you run GCP_FoldChange?')}
  if (!all(GCPlist$proteinINFO$FCcomparison[which(!is.na(GCPlist$proteinINFO$FCcomparison))] == GCPlist$proteinINFO$FCcomparison[which(!is.na(GCPlist$proteinINFO$FCcomparison))][1])) {stop('"FCcomparison" should have the same element in all rows!')}

  if (length(FDR)!=1) {stop("FDR must be exclusively TRUE or FALSE")}
  if (!is.logical(FDR)) {stop("FDR must be exclusively TRUE or FALSE")}
  if (is.na(FDR)) {stop("FDR must be exclusively TRUE or FALSE")}

  if (FDR) {
    if (!"ttest_PvaluesFDR" %in% colnames(GCPlist$proteinINFO)) {stop('"ttest_PvaluesFDR" is not present among the columns of proteinINFO. Did you run GCP_ttest?')}
  } else {
    if (!"ttest_Pvalues" %in% colnames(GCPlist$proteinINFO)) {stop('"ttest_Pvalues" is not present among the columns of proteinINFO. Did you run GCP_ttest?')}
  }


  if (length(log_base)!=1) {stop("log_base must be a number")}
  if (is.na(log_base)) {stop("log_base must be a number")}
  if (!is.numeric(log_base)) {stop("log_base must be a number")}

  if (!is.null(pcutoff_colored) & is.null(name_column_proteingroups)) {
    if (length(pcutoff_colored)!=1) {stop("pcutoff_colored must be a single number between 0 and 1")}
    if (is.na(pcutoff_colored)) {stop("pcutoff_colored must be a single number between 0 and 1, and not a missing value")}
    if (!is.numeric(pcutoff_colored)) {stop("pcutoff_colored must be a single number between 0 and 1")}
    if (pcutoff_colored>1 | pcutoff_colored<0) {stop("pcutoff_colored must be a single number between 0 and 1")}
  }

  if (!is.null(pcutoff_line)) {
    if (length(pcutoff_line)!=1) {stop("pcutoff_line must be a single number between 0 and 1")}
    if (is.na(pcutoff_line)) {stop("pcutoff_line must be a single number between 0 and 1, and not a missing value")}
    if (!is.numeric(pcutoff_line)) {stop("pcutoff_line must be a single number between 0 and 1")}
    if (pcutoff_line>1 | pcutoff_line<0) {stop("pcutoff_line must be a single number between 0 and 1")}
  }

  if (!is.null(pcutoff_prot_label) & !is.null(name_column_proteinlabels)) {
    if (length(pcutoff_prot_label)!=1) {stop("pcutoff_prot_label must be a single number between 0 and 1")}
    if (is.na(pcutoff_prot_label)) {stop("pcutoff_prot_label must be a single number between 0 and 1, and not a missing value")}
    if (!is.numeric(pcutoff_prot_label)) {stop("pcutoff_prot_label must be a single number between 0 and 1")}
    if (pcutoff_prot_label>1 | pcutoff_prot_label<0) {stop("pcutoff_prot_label must be a single number between 0 and 1")}

    if (length(name_column_proteinlabels)!=1) {stop("name_column_proteinlabels must be a character of length 1")}
    if (!is.character(name_column_proteinlabels)) {stop("name_column_proteinlabels must be a character of length 1")}
    if (is.na(name_column_proteinlabels)) {stop("name_column_proteinlabels must be a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$proteinINFO) == name_column_proteinlabels)) != 1) {stop("The name passed in name_column_proteinlabels must be a name of a column of the proteinINFO dataframe")}
  }

  if (!is.null(name_column_proteingroups)) {
    if (length(name_column_proteingroups)!=1) {stop("name_column_proteingroups must be a character of length 1")}
    if (!is.character(name_column_proteingroups)) {stop("name_column_proteingroups must be a character of length 1")}
    if (is.na(name_column_proteingroups)) {stop("name_column_proteingroups must be a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$proteinINFO) == name_column_proteingroups)) != 1) {stop("The name passed in name_column_proteingroups must be a name of a column of the proteinINFO dataframe")}

    if (!is.factor(pull(GCPlist$proteinINFO, name_column_proteingroups))) {
      GCPlist$proteinINFO[, name_column_proteingroups] <- as.factor(pull(GCPlist$proteinINFO, name_column_proteingroups))
    }
  }

  are.colors <- function (vect) {
    map_lgl(vect, ~tryCatch({
      is.matrix(col2rgb(.)) & ncol(col2rgb(.))>0
    }, error = function(e) {
      FALSE
    }))
  }

  if (!is.null(pcutoff_colored) & is.null(name_column_proteingroups)) {
    if (!is.character(col_pal_difference)) {stop("col_pal_difference must be a character vector")}
    if (length(col_pal_difference)!=3) {stop("col_pal_difference must have exactly 3 colors")}
    if (any(is.na(col_pal_difference))) {stop("col_pal_difference must not contain NAs")}
    if (!all(are.colors(col_pal_difference))) {stop(paste0('col_pal_difference must contain valid colors. The following are not:\n "',
                                                           paste0(col_pal_difference[which(!are.colors(col_pal_difference))], collapse = '", "'), '"'))}
    if (is.null(names(col_pal_difference))) {
      names(col_pal_difference) <- c("not significant", "decreased", "increased")
    } else {
      if (!all(names(col_pal_difference)%in%c("not significant", "decreased", "increased")) | !all(c("not significant", "decreased", "increased")%in%names(col_pal_difference))) {
        stop('if present, names of col_pal_difference must be: "not significant", "decreased", "increased"')
      }
    }
  }

  if (!is.null(name_column_proteingroups)) {
    if (!is.null(col_pal_groups)) {
      if (!is.character(col_pal_groups)) {stop("col_pal_groups must be a character vector")}
      if (any(is.na(col_pal_groups))) {stop("col_pal_groups must not contain NAs")}
    } else {
      col_pal_groups <- build_long_vector_of_colors()
    }

    if (length(col_pal_groups)<length(levels(pull(GCPlist$proteinINFO, name_column_proteingroups)))) {
      stop(paste0("There are ", length(levels(pull(GCPlist$proteinINFO, name_column_proteingroups))), " groups, and you have specified only ", length(col_pal_groups), " colors in col_pal_groups"))
    } else {
      col_pal_groups <- col_pal_groups[1:length(levels(pull(GCPlist$proteinINFO, name_column_proteingroups)))]
    }

    if (!all(are.colors(col_pal_groups))) {stop(paste0('col_pal_groups must contain valid colors. The following are not:\n "',
                                                           paste0(col_pal_groups[which(!are.colors(col_pal_groups))], collapse = '", "'), '"'))}
    if (is.null(names(col_pal_groups))) {
      names(col_pal_groups) <- levels(pull(GCPlist$proteinINFO, name_column_proteingroups))
    } else {
      if (any(duplicated(names(col_pal_groups)))) {"col_pal_groups has some duplicated in the names"}
      if (!all(names(col_pal_groups) %in% levels(pull(GCPlist$proteinINFO, name_column_proteingroups)) & levels(pull(GCPlist$proteinINFO, name_column_proteingroups)) %in% names(col_pal_groups))) {stop("the names of col_pal_groups don't correspond to the levels of the indicated groups of protein")}
    }
  }




  Volcano_tibble <- tibble(protid = GCPlist$proteinINFO$protid,
                           logFC = GCPlist$proteinINFO$logFC)

  if (FDR) {
    Volcano_tibble <- mutate(Volcano_tibble, pval = GCPlist$proteinINFO$ttest_PvaluesFDR)
  } else {
    Volcano_tibble <- mutate(Volcano_tibble, pval = GCPlist$proteinINFO$ttest_Pvalues)
  }

  Volcano_tibble <- mutate(Volcano_tibble, minusLogpval = -log10(pval))

  Volcano_tibble <- filter(Volcano_tibble, !is.na(logFC) & !is.na(minusLogpval))


  Volcano_plot <- ggplot(data = Volcano_tibble, mapping = aes(x = logFC, y = minusLogpval)) +
    geom_point()

  if (!is.null(pcutoff_colored) & is.null(name_column_proteingroups)) {
    Volcano_tibble <- mutate(Volcano_tibble, Difference = factor(ifelse(pval>=pcutoff_colored, "not significant",
                                                                        ifelse(logFC<0, "decreased",
                                                                               ifelse(logFC>0, "increased", NA))),
                                                                 levels = c("not significant", "decreased", "increased")))

    DifferenceLabels <- c(paste0("not significant:\n ", sum(Volcano_tibble$Difference=="not significant")),
                          paste0("decreased:\n ", sum(Volcano_tibble$Difference=="decreased")),
                          paste0("increased:\n ", sum(Volcano_tibble$Difference=="increased")))
    names(DifferenceLabels) <- c("not significant", "decreased", "increased")


    Volcano_plot <- ggplot(data = Volcano_tibble, mapping = aes(x = logFC, y = minusLogpval, col = Difference)) +
      geom_point() +
      scale_color_manual(values = col_pal_difference,
                         labels = DifferenceLabels)
  }

  if (!is.null(name_column_proteingroups)) {
    Volcano_tibble <- left_join(Volcano_tibble, GCPlist$proteinINFO[, c("protid", name_column_proteingroups)], by = "protid")

    Volcano_plot <- ggplot(data = Volcano_tibble, mapping = aes(x = logFC, y = minusLogpval, col = !!sym(name_column_proteingroups))) +
      geom_point() +
      scale_color_manual(values = col_pal_groups)
  }


  if (!is.null(pcutoff_line)) {
    Volcano_plot <- Volcano_plot +
      geom_hline(yintercept = -log10(pcutoff_line), linetype="dashed")
  }

  if (!is.null(pcutoff_prot_label) & !is.null(name_column_proteinlabels)) {
    Volcano_tibble <- left_join(Volcano_tibble, GCPlist$proteinINFO[, c("protid", name_column_proteinlabels)], by = "protid")

    Volcano_plot <- Volcano_plot +
      geom_text_repel(data=Volcano_tibble[which(pull(Volcano_tibble, pval) < pcutoff_prot_label),],
                      aes(x = logFC, y = minusLogpval, label = !!sym(name_column_proteinlabels)))
  }


  Volcano_plot <- Volcano_plot +
    theme_bw() +
    ggtitle(paste0("Changes in ", GCPlist$proteinINFO$FCcomparison[which(!is.na(GCPlist$proteinINFO$FCcomparison))][1])) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5))

  if (FDR) {
    Volcano_plot <- Volcano_plot +
      ylab("-log10(t-test FDR p-values)")
  } else {
    Volcano_plot <- Volcano_plot +
      ylab("-log10(t-test p-values)")
  }

  if (log_base == exp(1)) {
    Volcano_plot <- Volcano_plot +
      xlab("ln(Fold Change)")
  } else {
    Volcano_plot <- Volcano_plot +
      xlab(paste0("log", log_base, "(Fold Change)"))
  }


  return(Volcano_plot)
}
