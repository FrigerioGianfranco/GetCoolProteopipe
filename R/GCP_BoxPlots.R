#' Get Box plots from intensities.
#'
#' It creates box plots showing the distribution of proteins in each group.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param by_samples logical. If TRUE, a single box plot graph will be returned with a box for each sample. If FALSE, multiple box plots can be generated for each protein, as specified in the only_these_protid or only_the_first arguments.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The box plots will be performed only in the specified data.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used. To see the colors, run build_long_vector_of_colors().
#' @param Title NULL or character of length 1. The title you want to ad on the top of the plot. Used only if by_samples is TRUE.
#' @param only_these_protid NULL or a character vector of protid. Box plots will be performed only for those proteins. Used only if by_samples is FALSE.
#' @param only_the_first NULL or a numeric integer. If only_these_protid is NULL, box plots will be performed only for those first proteins. Used only if by_samples is FALSE.
#'
#' @return A ggplot object, if by_samples is TRUE; or a list of ggplot objects, if by_samples is FALSE.
#'
#' @export
GCP_BoxPlots <- function(GCPlist, by_samples = TRUE, name_column_groups = NULL, raw_or_LFQ = c("raw", "LFQ"), col_pal = NULL, Title = "Distribution of intensities", only_these_protid = NULL, only_the_first = NULL) {

  checkGCPlist(GCPlist)

  if (!is.logical(by_samples)) {stop("by_samples must be either TRUE or FALSE")}
  if (length(by_samples) != 1) {stop("by_samples must be either TRUE or FALSE")}
  if (is.na(by_samples)) {stop("by_samples must be either TRUE or FALSE")}

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

  if (!identical(tolower(raw_or_LFQ), c("raw", "lfq"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ_or_both must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("raw", "lfq"))

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


  named_colors_bygroup <- col_pal
  if (is.null(names(named_colors_bygroup))) {
    names(named_colors_bygroup) <- levels(pull(GCPlist$sampleINFO, name_column_groups))
  } else {
    if (any(duplicated(names(named_colors_bygroup)))) {"col_pal has some duplicated in the names"}
    if (!all(names(named_colors_bygroup) %in% levels(pull(GCPlist$sampleINFO, name_column_groups)) & levels(pull(GCPlist$sampleINFO, name_column_groups)) %in% names(named_colors_bygroup))) {stop("the names of col_pal don't correspond to the levels of name_column_groups")}
  }



  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw)
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ)
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }

  df_intensities_wg <- add_column(df_intensities,
                                  allwiththis = factor(NA, levels = levels(pull(GCPlist$sampleINFO, name_column_groups))),
                                  .after = 1)
  colnames(df_intensities_wg)[2] <- name_column_groups

  for (i in 1:length(pull(df_intensities_wg, 1))) {
    df_intensities_wg[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
  }

  if (by_samples) {

    if (!is.null(Title)) {
      if (length(Title)!=1) {
        stop("Pass a suitable Title!")
      } else {
        if (!is.character(Title) & !is.numeric(Title) & !is.logical(Title)) {
          stop("Pass a suitable Title!")
        }
      }
    }


    df_intensities_wg_long <- pivot_longer(df_intensities_wg,
                                           cols = colnames(df_intensities_wg)[which(!colnames(df_intensities_wg) %in% c("samples", name_column_groups))],
                                           names_to = "protid",
                                           values_to = "intensity")

    df_intensities_wg_long$samples <- factor(df_intensities_wg_long$samples, levels = unique(df_intensities_wg_long$samples))


    the_plot <- ggplot(data = df_intensities_wg_long, aes(y = intensity, x = samples, fill = !!sym(name_column_groups))) +
      geom_boxplot(colour = "black") +
      scale_fill_manual(values = named_colors_bygroup) +
      ggtitle(Title) +
      theme_bw() +
      theme(plot.title = element_text(face="bold", hjust = 0.5))

    return(the_plot)

  } else {

    if (!is.null(only_these_protid)) {
      if (!is.character(only_these_protid)) {stop("if not NULL, only_these_protid must be a character vector")}
      if (any(is.na(only_these_protid))) {stop("if not NULL, only_these_protid must not contain NAs")}
      if (!all(only_these_protid %in% colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])) {stop("if not NULL, only_these_protid must contain names of columns of the data intensities")}
    }

    if (!is.null(only_the_first)) {
      if (length(only_the_first)!=1) {stop("if not NULL, only_the_first must be a numeric of lenght 1")}
      if (!is.numeric(only_the_first)) {stop("if not NULL, only_the_first must be a numeric of lenght 1")}
      if (only_the_first<1) {stop("if not NULL, only_the_first must be a numeric integer greater or equal to 1")}
      if (only_the_first != as.integer(only_the_first)) {stop("if not NULL, only_the_first must be a numeric integer greater or equal to 1")}
      if (is.na(only_the_first)) {stop("if not NULL, only_the_first must be a numeric of lenght 1, not a missing value")}
    }



    if (is.null(only_these_protid) & is.null(only_the_first)) {
      prot_to_use <- colnames(df_intensities)[-which(colnames(df_intensities) == "samples")]
    } else if (!is.null(only_these_protid)) {
      prot_to_use <- only_these_protid
    } else if (!is.null(only_the_first)) {
      if (only_the_first > length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])) {
        only_the_first <- length(colnames(df_intensities)[-which(colnames(df_intensities) == "samples")])
      }
      prot_to_use <- colnames(df_intensities)[-which(colnames(df_intensities) == "samples")][1:only_the_first]
    }



    boxplots_list <- vector(mode = "list", length = length(prot_to_use))
    names(boxplots_list) <- prot_to_use


    for (this_prot in prot_to_use) {

      index_in_the_INFO <- which(GCPlist$proteinINFO$protid == this_prot)



      this_plot <- ggplot(data = df_intensities_wg, aes(y = !!sym(this_prot), x = !!sym(name_column_groups), fill = !!sym(name_column_groups))) +
        geom_boxplot(colour = "black") +
        scale_fill_manual(values = named_colors_bygroup) +
        theme_bw()


      if ("Protein names" %in% colnames(GCPlist$proteinINFO)) {
        if (!is.na(GCPlist$proteinINFO$`Protein names`[index_in_the_INFO])) {

          this_plot <- this_plot +
            labs(title = GCPlist$proteinINFO$`Protein names`[index_in_the_INFO]) +
            theme(plot.title = element_text(hjust = 0.5))
        }
      }

      if ("Accession" %in% colnames(GCPlist$proteinINFO)) {
        if (!is.na(GCPlist$proteinINFO$Accession[index_in_the_INFO])) {

          this_plot <- this_plot +
            labs(subtitle = GCPlist$proteinINFO$Accession[index_in_the_INFO]) +
            theme(plot.subtitle = element_text(hjust = 0.5))
        }
      }

      boxplots_list[[which(prot_to_use == this_prot)]] <- this_plot

    }

    return(boxplots_list)

  }
}


