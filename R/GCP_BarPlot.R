#' Get bar plots with protein numbers.
#'
#' It create bar plots to show the number of valid proteins per samples.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param raw_or_LFQ_or_both one of the following: "raw", "LFQ", "both". The barplot will be performed only in the specified data, or a facet_wrap plot will be generated if "both" is specified.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used (see function build_long_vector_of_colors).
#' @param label_numbers logical. If TRUE, it adds the number of proteins on the top of each bar.
#' @param showCV logical. If TRUE, it adds to the graphs the coefficient of variation (CV%) of the number of proteins for each group.
#' @param rotate_sample_names logical. If TRUE, it rotate sample name labels of x-axis in vertical position.
#'
#' @return A ggplot object.
#'
#' @export
GCP_BarPlot <- function(GCPlist, name_column_groups = NULL, raw_or_LFQ_or_both = c("raw", "LFQ", "both"), col_pal = NULL, label_numbers = TRUE, showCV = FALSE, rotate_sample_names = FALSE) {

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

  if (!identical(tolower(raw_or_LFQ_or_both), c("raw", "lfq", "both"))) {
    if (length(raw_or_LFQ_or_both) != 1) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ", "both"')}
    if (is.na(raw_or_LFQ_or_both)) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ", "both"')}
  }
  raw_or_LFQ_or_both <- tolower(raw_or_LFQ_or_both)
  raw_or_LFQ_or_both <- match.arg(raw_or_LFQ_or_both, c("raw", "lfq", "both"))

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

  if (!is.logical(label_numbers)) {stop("label_numbers must be either TRUE or FALSE")}
  if (length(label_numbers) != 1) {stop("label_numbers must be either TRUE or FALSE")}
  if (is.na(label_numbers)) {stop("label_numbers must be either TRUE or FALSE")}

  if (!is.logical(showCV)) {stop("showCV must be either TRUE or FALSE")}
  if (length(showCV) != 1) {stop("showCV must be either TRUE or FALSE")}
  if (is.na(showCV)) {stop("showCV must be either TRUE or FALSE")}

  if (!is.logical(rotate_sample_names)) {stop("rotate_sample_names must be either TRUE or FALSE")}
  if (length(rotate_sample_names) != 1) {stop("rotate_sample_names must be either TRUE or FALSE")}
  if (is.na(rotate_sample_names)) {stop("rotate_sample_names must be either TRUE or FALSE")}


  if (raw_or_LFQ_or_both == "raw") {

    prot_num_df_raw <- summarise_at(GCPlist$quant_raw, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_raw_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "proteins"), prot_num_df_raw), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_raw <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_raw_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_raw[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_raw, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_raw, colnames(GCPlist$sampleINFO)[1])))

    the_barplot <- ggplot(data = df_barplot_raw, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["proteins"]], fill = .data[[name_column_groups]])) +
      geom_col() +
      ggtitle("Protein numbers in quant_raw") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    if (showCV) {
      summary_for_CV_raw <- df_barplot_raw %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(proteins), SD = sd(proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_raw <- summary_for_CV_raw$CV
      names(CV_raw) <- pull(summary_for_CV_raw, name_column_groups)

      CV_labels <- paste0(names(CV_raw), "\n CV: ",  CV_raw, "%")
      names(CV_labels) <- names(CV_raw)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" CV: ",  CV_raw, "%")
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
      max_proteins <- max(df_barplot_raw$proteins)
      the_barplot <- the_barplot + geom_text(aes(label = .data[["proteins"]], y = .data[["proteins"]]+max_proteins*0.022))
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}

  } else if (raw_or_LFQ_or_both == "lfq") {

    prot_num_df_LFQ <- summarise_at(GCPlist$quant_LFQ, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_LFQ_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "proteins"), prot_num_df_LFQ), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_LFQ <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_LFQ_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_LFQ[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_LFQ, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_LFQ, colnames(GCPlist$sampleINFO)[1])))

    the_barplot <- ggplot(data = df_barplot_LFQ, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["proteins"]], fill = .data[[name_column_groups]])) +
      geom_col() +
      ggtitle("Protein numbers in quant_LFQ") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    if (showCV) {
      summary_for_CV_LFQ <- df_barplot_LFQ %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(proteins), SD = sd(proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_LFQ <- summary_for_CV_LFQ$CV
      names(CV_LFQ) <- pull(summary_for_CV_LFQ, name_column_groups)

      CV_labels <- paste0(names(CV_LFQ), "\n CV: ",  CV_LFQ, "%")
      names(CV_labels) <- names(CV_LFQ)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" CV: ",  CV_LFQ, "%")
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
      max_proteins <- max(df_barplot_LFQ$proteins)
      the_barplot <- the_barplot + geom_text(aes(label = .data[["proteins"]], y = .data[["proteins"]]+max_proteins*0.022))
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}

  } else if (raw_or_LFQ_or_both == "both") {

    prot_num_df_raw <- summarise_at(GCPlist$quant_raw, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_raw_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "proteins"), prot_num_df_raw), name_first_column = colnames(GCPlist$sampleINFO)[1])

    prot_num_df_LFQ <- summarise_at(GCPlist$quant_LFQ, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!="protid")], ~ sum(!is.na(.x)))

    prot_num_df_LFQ_t <- GetFeatistics::transpose_feat_table(bind_cols(tibble(colcont = "proteins"), prot_num_df_LFQ), name_first_column = colnames(GCPlist$sampleINFO)[1])

    df_barplot_raw <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_raw_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_LFQ <- left_join(x = GCPlist$sampleINFO, y = prot_num_df_LFQ_t, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_prot_summ"))

    df_barplot_both <- bind_rows(list(raw = df_barplot_raw,
                                      LFQ = df_barplot_LFQ),
                                 .id = "intensity_table")
    df_barplot_both$intensity_table <- factor(df_barplot_both$intensity_table, levels = c("raw", "LFQ"))

    df_barplot_both[,colnames(GCPlist$sampleINFO)[1]] <- factor(pull(df_barplot_both, colnames(GCPlist$sampleINFO)[1]), levels = unique(pull(df_barplot_both, colnames(GCPlist$sampleINFO)[1])))

    the_barplot <- ggplot(data = df_barplot_both, aes(x = .data[[colnames(GCPlist$sampleINFO)[1]]], y = .data[["proteins"]], fill = .data[[name_column_groups]])) +
      geom_col() +
      facet_wrap(~ intensity_table) +
      ggtitle("Protein numbers") +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5))

    if (showCV) {
      summary_for_CV_raw <- df_barplot_raw %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(proteins), SD = sd(proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_raw <- summary_for_CV_raw$CV
      names(CV_raw) <- pull(summary_for_CV_raw, name_column_groups)

      summary_for_CV_LFQ <- df_barplot_LFQ %>%
        group_by(!!sym(name_column_groups)) %>%
        summarise(Mean = mean(proteins), SD = sd(proteins)) %>%
        mutate(CV = round((SD/Mean)*100, digits = 1))

      CV_LFQ <- summary_for_CV_LFQ$CV
      names(CV_LFQ) <- pull(summary_for_CV_LFQ, name_column_groups)


      CV_labels <- paste0(names(CV_raw), "\n CV raw: ",  CV_raw, "%\n CV LFQ: ", CV_LFQ, "%")
      names(CV_labels) <- names(CV_raw)

      if (name_column_groups == "allwiththis") {
        CV_labels <- paste0(" CV raw: ",  CV_raw, "%\n CV LFQ: ", CV_LFQ, "%")
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
      max_proteins <- max(df_barplot_both$proteins)
      the_barplot <- the_barplot + geom_text(aes(label = .data[["proteins"]], y = .data[["proteins"]]+max_proteins*0.022))
    }

    if (name_column_groups == "allwiththis" & !showCV) {the_barplot <- the_barplot + theme(legend.position = "none")}


  } else {
      stop('raw_or_LFQ must be one of "raw", "LFQ", "both"')
  }



  if (rotate_sample_names) {
    the_barplot <- the_barplot +
      theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust = 0))
  }

  return(the_barplot)
}



