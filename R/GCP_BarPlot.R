#' Get bar plots with protein numbers.
#'
#' It create bar plots to show the number of valid proteins per samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param raw_or_LFQ_or_both one of the following: "raw", "LFQ", "both". The barplot will be performed only in the specified data, or a facet_wrap plot will be generated if "both" is specified.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used (see function build_long_vector_of_colors).
#' @param bar_width NULL or a number. You can pass here the width of the bars.
#' @param label_numbers logical. If TRUE, it adds the number of proteins on the top of each bar.
#' @param label_numbers_size NULL or numeric of length 1. If specified and if label_numbers is TRUE, this is the size of the numbers of proteins of the top of each bar.
#' @param showCV logical. If TRUE, it adds to the graphs the coefficient of variation (CV%) of the number of proteins for each group.
#' @param rotate_sample_names logical. If TRUE, it rotate sample name labels of x-axis in vertical position.
#'
#' @return A ggplot object.
#'
#' @export
GCP_BarPlot <- function(GCPlist, name_column_groups = NULL, raw_or_LFQ_or_both = getOption("GetCoolProteopipe.raw_or_LFQ"), col_pal = NULL, bar_width = NULL, label_numbers = TRUE, label_numbers_size = NULL, showCV = FALSE, rotate_sample_names = FALSE) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

  } else {
    if ("allwiththis" %in% colnames(GCPlist$sampleINFO)) {stop("Please, don't call a column of the sampleINFO data frame 'allwiththis' as I need to create one with this name now")}
    GCPlist$sampleINFO <- mutate(GCPlist$sampleINFO, allwiththis = as.factor("theOnlyGroup"))

    name_column_groups <- "allwiththis"
  }

  if (!identical(tolower(raw_or_LFQ_or_both), c("lfq", "raw", "both"))) {
    if (length(raw_or_LFQ_or_both) != 1) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ", "both"')}
    if (is.na(raw_or_LFQ_or_both)) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ", "both"')}
  }
  raw_or_LFQ_or_both <- tolower(raw_or_LFQ_or_both)
  raw_or_LFQ_or_both <- match.arg(raw_or_LFQ_or_both, c("lfq", "raw", "both"))

  if (raw_or_LFQ_or_both == "lfq") {
    cat("\n -- LFQ data are used --\n\n")
  } else if (raw_or_LFQ_or_both == "raw") {
    cat("\n -- raw data are used --\n\n")
  } else if (raw_or_LFQ_or_both == "both") {
    cat("\n -- both raw and LFQ data are compared --\n\n")
  }


  if (!is.null(col_pal)) {
    if (!is.character(col_pal)) stop("col_pal must be a character vector")
    if (any(is.na(col_pal))) stop("col_pal must not contain NAs")
  } else {
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


  if (raw_or_LFQ_or_both == "raw") {

    prot_num_df_raw <- summarise_at(GCPlist$quant_raw, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_raw_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "Proteins"), prot_num_df_raw), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_raw <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_raw_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_raw[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_raw, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_raw, colnames(GCPlist$sampleINFO)[1])))


    if (is.null(bar_width)) {
      the_barplot <- ggplot(data = df_barplot_raw, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col() +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      the_barplot <- ggplot(data = df_barplot_raw, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col(width = bar_width) +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    }



    if (showCV) {
      summary_for_CV_raw <- df_barplot_raw %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(Proteins), SD = sd(Proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_raw <- summary_for_CV_raw$CV
      names(CV_raw) <- pull(summary_for_CV_raw, name_column_groups)

      CV_labels <- paste0(names(CV_raw), " (CV: ",  CV_raw, "%)")
      names(CV_labels) <- names(CV_raw)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" (CV: ",  CV_raw, "%)")
        names(CV_labels) <- names(CV_raw)
      }

      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups,
                          labels = CV_labels)

    } else {
      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups)
    }

    if (label_numbers) {
      max_proteins <- max(df_barplot_raw$Proteins)

      if (is.null(label_numbers_size)) {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022))
      } else {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022), size = label_numbers_size)
      }
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}

  } else if (raw_or_LFQ_or_both == "lfq") {

    prot_num_df_LFQ <- summarise_at(GCPlist$quant_LFQ, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_LFQ_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "Proteins"), prot_num_df_LFQ), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_LFQ <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_LFQ_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_LFQ[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_LFQ, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_LFQ, colnames(GCPlist$sampleINFO)[1])))

    if (is.null(bar_width)) {
      the_barplot <- ggplot(data = df_barplot_LFQ, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col() +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      the_barplot <- ggplot(data = df_barplot_LFQ, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col(width = bar_width) +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    }


    if (showCV) {
      summary_for_CV_LFQ <- df_barplot_LFQ %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(Proteins), SD = sd(Proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_LFQ <- summary_for_CV_LFQ$CV
      names(CV_LFQ) <- pull(summary_for_CV_LFQ, name_column_groups)

      CV_labels <- paste0(names(CV_LFQ), " (CV: ",  CV_LFQ, "%)")
      names(CV_labels) <- names(CV_LFQ)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" (CV: ",  CV_LFQ, "%)")
        names(CV_labels) <- names(CV_LFQ)
      }

      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups,
                          labels = CV_labels)

    } else {
      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups)
    }


    if (label_numbers) {
      max_proteins <- max(df_barplot_LFQ$Proteins)

      if (is.null(label_numbers_size)) {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022))
      } else {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022), size = label_numbers_size)
      }
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}

  } else if (raw_or_LFQ_or_both == "both") {

    prot_num_df_raw <- summarise_at(GCPlist$quant_raw, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_raw_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "Proteins"), prot_num_df_raw), name_first_column = colnames(GCPlist$sampleINFO)[1])

    prot_num_df_LFQ <- summarise_at(GCPlist$quant_LFQ, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_LFQ_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "Proteins"), prot_num_df_LFQ), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_raw <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_raw_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_LFQ <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_LFQ_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_both <- bind_rows(list(raw = df_barplot_raw,
                                      LFQ = df_barplot_LFQ),
                                 .id = "intensity_table")
    df_barplot_both$intensity_table <- factor(df_barplot_both$intensity_table, levels = c("raw", "LFQ"))

    df_barplot_both[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_both, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_both, colnames(GCPlist$sampleINFO)[1])))

    if (is.null(bar_width)) {
      the_barplot <- ggplot(data = df_barplot_both, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col() +
        facet_wrap(~ intensity_table) +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    } else {
      the_barplot <- ggplot(data = df_barplot_both, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["Proteins"]], fill = .data[[name_column_groups]])) +
        geom_col(width = bar_width) +
        facet_wrap(~ intensity_table) +
        ggtitle("Number of proteins per sample") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5))
    }

    if (showCV) {
      summary_for_CV_raw <- df_barplot_raw %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(Proteins), SD = sd(Proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_raw <- summary_for_CV_raw$CV
      names(CV_raw) <- pull(summary_for_CV_raw, name_column_groups)

      summary_for_CV_LFQ <- df_barplot_LFQ %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(Proteins), SD = sd(Proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_LFQ <- summary_for_CV_LFQ$CV
      names(CV_LFQ) <- pull(summary_for_CV_LFQ, name_column_groups)


      CV_labels <- paste0(names(CV_raw), " (CV raw: ",  CV_raw, "%) (CV LFQ: ", CV_LFQ, "%)")
      names(CV_labels) <- names(CV_raw)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" (CV raw: ",  CV_raw, "%) (CV LFQ: ", CV_LFQ, "%)")
        names(CV_labels) <- names(CV_raw)
      }


      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups,
                          labels = CV_labels)

    } else {
      the_barplot <- the_barplot +
        scale_fill_manual(values = colors_of_groups)
    }


    if (label_numbers) {
      max_proteins <- max(df_barplot_both$Proteins)

      if (is.null(label_numbers_size)) {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022))
      } else {
        the_barplot <- the_barplot + geom_text(aes(label = .data[["Proteins"]], y = .data[["Proteins"]]+max_proteins*0.022), size = label_numbers_size)
      }
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}


  } else {
      stop('raw_or_LFQ must be one of "raw", "LFQ", "both"')
  }



  if (rotate_sample_names) {
    the_barplot <- the_barplot +
      theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
  }


  if (raw_or_LFQ_or_both == "raw") {

    if (any(map_lgl(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")], ~ any(isTRUE(.x == 0))))) {
      warning("There are some zeros in the raw table, are you sure you don't want to apply before the function GCP_ReplaceZerowithNA ?")
    } else if (!any(map_lgl(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")], ~ any(is.na(.x))))) {
      warning("There are literally no NAs in the entire raw table, that's why all the numbers are the same!")
    }



  } else if (raw_or_LFQ_or_both == "lfq") {
    if (any(map_lgl(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")], ~ any(isTRUE(.x == 0))))) {
      warning("There are some zeros in the LFQ table, are you sure you don't want to apply before the function GCP_ReplaceZerowithNA ?")
    } else if (!any(map_lgl(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")], ~ any(is.na(.x))))) {
      warning("There are literally no NAs in the entire LFQ table, that's why all the numbers are the same!")
    }


  } else if (raw_or_LFQ_or_both == "both") {
    if (any(c(any(map_lgl(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")], ~ any(isTRUE(.x == 0)))), any(map_lgl(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")], ~ any(isTRUE(.x == 0))))))) {
      warning("There are some zeros in the intensity tables, are you sure you don't want to apply before the function GCP_ReplaceZerowithNA ?")
    } else if (!any(c(any(map_lgl(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")], ~ any(is.na(.x)))), any(map_lgl(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")], ~ any(is.na(.x))))))) {
      warning("There are literally no NAs in the intensity tables, that's why all the numbers are the same!")
    }
  }


  return(the_barplot)
}



