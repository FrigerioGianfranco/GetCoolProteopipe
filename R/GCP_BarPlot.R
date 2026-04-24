#' Get bar plots with protein numbers.
#'
#' It create bar plots to show the number of valid proteins per samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used (see function build_long_vector_of_colors).
#' @param bar_width NULL or a number. You can pass here the width of the bars.
#' @param label_numbers logical. If TRUE, it adds the number of proteins on the top of each bar.
#' @param label_numbers_size NULL or numeric of length 1. If specified and if label_numbers is TRUE, this is the size of the numbers of proteins of the top of each bar.
#' @param showCV logical. If TRUE, it adds to the graphs the coefficient of variation (CV%) of the number of proteins for each group.
#' @param rotate_sample_names logical. If TRUE, it rotate sample name labels of x-axis in vertical position.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#' \dontrun{
#'
#' Fig07_BarPlot_post_filteringNA <- GCP_BarPlot(GCPlist = GCPlist07,
#'                                               label_numbers = TRUE,
#'                                               showCV = TRUE,
#'                                               rotate_sample_names = TRUE)
#' export_figures(Fig07_BarPlot_post_filteringNA)
#'
#' }
#'
#'
#'
#' @export
GCP_BarPlot <- function(GCPlist, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), col_pal = getOption("GetCoolProteopipe.col_pal"), bar_width = NULL, label_numbers = TRUE, label_numbers_size = NULL, showCV = FALSE, rotate_sample_names = FALSE) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }
    if (any(table(pull(GCPlist$sampleINFO, name_column_groups)) == 0)) {
      GCPlist$sampleINFO[,name_column_groups] <- droplevels(pull(GCPlist$sampleINFO, name_column_groups))
    }

  } else {
    cat("\n -- The name_column_groups considered is NULL --\n\n")
    if ("allwiththis" %in% colnames(GCPlist$sampleINFO)) {stop("Please, don't call a column of the sampleINFO data frame 'allwiththis' as I need to create one with this name now")}
    GCPlist$sampleINFO <- mutate(GCPlist$sampleINFO, allwiththis = as.factor("theOnlyGroup"))

    name_column_groups <- "allwiththis"
  }

  if (!is.null(col_pal)) {
    if (!is.character(col_pal)) stop("col_pal must be a character vector")
    if (any(is.na(col_pal))) stop("col_pal must not contain NAs")
    if (is.null(names(col_pal))) {
      cat(paste0(' --- col_pal is  c("', paste0(col_pal, collapse = '", "'), '") ---\n\n'))
    } else {
      cat(' --- col_pal is  c(')
      for (i in seq(length(col_pal))) {
        if (names(col_pal)[i]!="") {
          cat(paste0(names(col_pal)[i], ' = "'))
        } else {
          cat('"')
        }
        cat(paste0(col_pal[i], '"'))
        if (i != length(col_pal)) {
          cat(', ')
        }
      }
      cat(') ---\n\n')
    }
  } else {
    cat(" -- col_pal is  NULL --\n\n")
    col_pal <- build_long_vector_of_colors()
  }

  if (length(col_pal)<length(levels(pull(GCPlist$sampleINFO, name_column_groups)))) {
    stop(paste0("There are ", length(levels(pull(GCPlist$sampleINFO, name_column_groups))), " groups, and you have specified only ", length(col_pal), " colors in col_pal"))
  } else {
    col_pal <- col_pal[1:length(levels(pull(GCPlist$sampleINFO, name_column_groups)))]
  }

  are.colors <- function (vect) {
    map_lgl(vect, ~tryCatch({
      is.matrix(col2rgb(.)) & ncol(col2rgb(.))>0
    }, error = function(e) {
      FALSE
    }))
  }
  if (!all(are.colors(col_pal))) {stop(paste0('col_pal must contain valid colors. In particular, these are not: "', paste(col_pal[!are.colors(col_pal)], collapse = '", "')), '"')}

  colors_of_groups <- col_pal
  if (is.null(names(colors_of_groups))) {
    names(colors_of_groups) <- levels(pull(GCPlist$sampleINFO, name_column_groups))
  } else {
    if (any(duplicated(names(colors_of_groups)))) {"col_pal has some duplicated in the names"}
    if (!all(names(colors_of_groups) %in% levels(pull(GCPlist$sampleINFO, name_column_groups)) & levels(pull(GCPlist$sampleINFO, name_column_groups)) %in% names(colors_of_groups))) {stop("the names of col_pal don't correspond to the levels of name_column_groups")}
  }

  if (!is.null(bar_width)) {
    if (length(bar_width)!=1) {stop("bar_width must be a numeric of length 1 (or NULL)")}
    if (!is.numeric(bar_width)) {stop("bar_width must be a numeric of length 1 (or NULL)")}
    if (is.na(bar_width)) {stop("bar_width must be a numeric of length 1 (or NULL), not a missing value")}
  }

  if (!is.logical(label_numbers)) {stop("label_numbers must be either TRUE or FALSE")}
  if (length(label_numbers) != 1) {stop("label_numbers must be either TRUE or FALSE")}
  if (is.na(label_numbers)) {stop("label_numbers must be either TRUE or FALSE")}

  if (label_numbers) {
    if (!is.null(label_numbers_size)) {
      if (length(label_numbers_size)!=1) {stop("label_numbers_size must be a numeric of length 1")}
      if (!is.numeric(label_numbers_size)) {stop("label_numbers_size must be a numeric of length 1")}
      if (is.na(label_numbers_size)) {stop("label_numbers_size must be a numeric of length 1, not a missing value")}
    }
  }

  if (!is.logical(showCV)) {stop("showCV must be either TRUE or FALSE")}
  if (length(showCV) != 1) {stop("showCV must be either TRUE or FALSE")}
  if (is.na(showCV)) {stop("showCV must be either TRUE or FALSE")}

  if (!is.logical(rotate_sample_names)) {stop("rotate_sample_names must be either TRUE or FALSE")}
  if (length(rotate_sample_names) != 1) {stop("rotate_sample_names must be either TRUE or FALSE")}
  if (is.na(rotate_sample_names)) {stop("rotate_sample_names must be either TRUE or FALSE")}

  prot_num_df <- summarise_at(GCPlist$intensities, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!="protid")], ~ sum(!is.na(.x)))

  prot_num_df_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "Proteins"), prot_num_df), name_first_column = colnames(GCPlist$sampleINFO)[1])

  df_barplot <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

  df_barplot[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot, colnames(GCPlist$sampleINFO)[1])))


  if (is.null(bar_width)) {
    the_barplot <- ggplot(data = df_barplot, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
      geom_col() +
      ggtitle("Number of proteins per sample") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    the_barplot <- ggplot(data = df_barplot, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
      geom_col(width = bar_width) +
      ggtitle("Number of proteins per sample") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))
  }

  if (showCV) {
    summary_for_CV <- df_barplot %>%
      group_by(!!sym(name_column_groups)) %>%
      summarise(Mean = mean(Proteins), SD = sd(Proteins)) %>%
      mutate(CV = round((SD/Mean)*100, digits = 1))

    CV <- summary_for_CV$CV
    names(CV) <- pull(summary_for_CV, name_column_groups)

    CV_labels <- paste0(names(CV), " (CV: ",  CV, "%)")
    names(CV_labels) <- names(CV)

    if (name_column_groups == "allwiththis") {
      CV_labels <- paste0(" (CV: ",  CV, "%)")
      names(CV_labels) <- names(CV)
    }

    the_barplot <- the_barplot +
      scale_fill_manual(values = colors_of_groups,
                        labels = CV_labels)

  } else {
    the_barplot <- the_barplot +
      scale_fill_manual(values = colors_of_groups)
  }

  if (label_numbers) {
    max_proteins <- max(df_barplot$Proteins)

    if (is.null(label_numbers_size)) {
      the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022))
    } else {
      the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022), size = label_numbers_size)
    }
  }

  if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}

  if (rotate_sample_names) {
    the_barplot <- the_barplot +
      theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
  }


  if (any(map_lgl(GCPlist$intensities[, which(colnames(GCPlist$intensities)!="protid")], ~ any(isTRUE(.x == 0))))) {
    warning("There are some zeros in the intensities table, are you sure you don't want to apply before the function GCP_ReplaceZerowithNA() ?")
  } else if (!any(map_lgl(GCPlist$intensities[, which(colnames(GCPlist$intensities)!="protid")], ~ any(is.na(.x))))) {
    warning("There are literally no NAs in the intensities table, that's why all the numbers are the same!")
  }

  return(the_barplot)
}



