#' Plot Principal Component analysis.
#'
#' It plots a principal component analysis.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The principal component analysis will be performed only in the specified data.
#' @param scores_or_loadings one of the following: "scores", "loadings". Specify here if you want to plot the scores or the loadings.
#' @param PC_to_plot character of length 2. Specify here the two principal components to plot.
#' @param center logical. Whether the variables should be shifted to be zero centered (as in the prcomp function).
#' @param scale. logical. whether the variables should be scaled to have unit variance before the analysis takes place (as in prcomp function).
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups. Specify it only if you want to color the points of the score plot.
#' @param name_column_labels NULL or character of length 1. The name of the column of the sampleINFO table containing the sample names. Specify it only if you want to add a label to the points of the score plot.
#' @param col_pal a character vector containing colors for groups for the score plot. If NULL, colors from the pals package will be used (see function GetFeatistics::build_long_vector_of_colors).
#' @param ellipses_on_score logical. If you specified name_column_groups and this is TRUE, ellipses will be added to the score plot.
#' @param name_column_groups_loading NULL or character of length 1. The name of the column of the proteinINFO table containing the protein groups. Specify it only if you want to color the points of the loading plot.
#' @param name_column_labels_loading NULL or character of length 1. The name of the column of the proteinINFO table containing the protein/protein names. Specify it only if you want to add a label to the points of the loading plot.
#' @param col_pal_loading a character vector containing colors for groups for the loading plot. If NULL, colors from the pals package will be used (see function GetFeatistics::build_long_vector_of_colors).
#' @param ellipses_on_loading logical. If you specified name_column_groups_loading and this is TRUE, ellipses will be added to the loading plot.
#'
#' @return A ggplot object.
#'
#' @import ggdendro
#'
#' @export
GCP_plotPCA <- function(GCPlist, raw_or_LFQ = getOption("GetCoolProteopipe.raw_or_LFQ"),
                        scores_or_loadings = c("scores", "loadings"), PC_to_plot = c("PC1", "PC2"),
                        center = TRUE, scale. = FALSE,
                        name_column_groups = NULL, name_column_labels = NULL, col_pal = NULL, ellipses_on_score = TRUE,
                        name_column_groups_loading = NULL, name_column_labels_loading = NULL, col_pal_loading = NULL, ellipses_on_loading = FALSE) {

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

  if (!identical(tolower(scores_or_loadings), c("scores", "loadings"))) {
    if (length(scores_or_loadings) != 1) {stop('scores_or_loadings must be one of "scores", "loadings"')}
    if (is.na(scores_or_loadings)) {stop('scores_or_loadings must be one of "scores", "loadings"')}
  }
  scores_or_loadings <- tolower(scores_or_loadings)
  scores_or_loadings <- match.arg(scores_or_loadings, c("scores", "loadings"))

  if (length(PC_to_plot) != 2) {stop("PC_to_plot must contain exactly 2 elements")}
  if (!is.character(PC_to_plot)) {stop("PC_to_plot must be a character vector")}
  if (any(is.na(PC_to_plot))) {stop("PC_to_plot must not contain NAs")}
  if (raw_or_LFQ == "raw") {
    if (PC_to_plot[1] %in% pull(GCPlist$quant_raw, 1)) {stop(paste0('This is a problem..: quant_raw already contains a protein named "', PC_to_plot[1], '"'))}
    if (PC_to_plot[2] %in% pull(GCPlist$quant_raw, 1)) {stop(paste0('This is a problem..: quant_raw already contains a protein named "', PC_to_plot[2], '"'))}
  } else if (raw_or_LFQ == "lfq") {
    if (PC_to_plot[1] %in% pull(GCPlist$quant_LFQ, 1)) {stop(paste0('This is a problem..: quant_LFQ already contains a protein named "', PC_to_plot[1], '"'))}
    if (PC_to_plot[2] %in% pull(GCPlist$quant_LFQ, 1)) {stop(paste0('This is a problem..: quant_LFQ already contains a protein named "', PC_to_plot[2], '"'))}
  }

  are.colors <- function (vect) {
    map_lgl(vect, ~tryCatch({
      is.matrix(col2rgb(.)) & ncol(col2rgb(.))>0
    }, error = function(e) {
      FALSE
    }))
  }
  used_the_long_vector_of_colors <- FALSE


  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}
    if (name_column_groups == PC_to_plot[1]) {stop(paste0('This is a problem..: name_column_groups should not be named as "', PC_to_plot[1], '"'))}
    if (name_column_groups == PC_to_plot[2]) {stop(paste0('This is a problem..: name_column_groups should not be named as "', PC_to_plot[2], '"'))}
    if (name_column_groups == "allwiththis") {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}
    if (name_column_groups == "thesearethesamplenamesused") {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_groups, thanks!")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

    if (!is.null(col_pal)) {
      if (!is.character(col_pal)) {stop("col_pal must be a character vector")}
      if (any(is.na(col_pal))) {stop("col_pal must not contain NAs")}

    } else {
      col_pal <- build_long_vector_of_colors()
      used_the_long_vector_of_colors <- TRUE
    }

    if (length(col_pal)<length(levels(pull(GCPlist$sampleINFO, name_column_groups)))) {
      stop(paste0("There are ", length(levels(pull(GCPlist$sampleINFO, name_column_groups))), " groups, and you have specified only ", length(col_pal), " colors in col_pal"))
    } else {
      col_pal <- col_pal[1:length(levels(pull(GCPlist$sampleINFO, name_column_groups)))]
    }

    if (!all(are.colors(col_pal))) {stop(paste0('col_pal must contain valid colors. In particular, these are not: "', paste(col_pal[!are.colors(col_pal)], collapse = '", "')), '"')}
  }

  if (length(ellipses_on_score)!=1) {stop("ellipses_on_score must be exclusively TRUE or FALSE")}
  if (!is.logical(ellipses_on_score)) {stop("ellipses_on_score must be exclusively TRUE or FALSE")}
  if (is.na(ellipses_on_score)) {stop("ellipses_on_score must be exclusively TRUE or FALSE")}


  if (!is.null(name_column_groups_loading)) {
    if (length(name_column_groups_loading)!=1) {stop("name_column_groups_loading must be NULL or a character of length 1")}
    if (!is.character(name_column_groups_loading)) {stop("name_column_groups_loading must be NULL or a character of length 1")}
    if (is.na(name_column_groups_loading)) {stop("name_column_groups_loading must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$proteinINFO) == name_column_groups_loading)) != 1) {stop("The name passed in name_column_groups_loading must be a name of a column of the proteinINFO dataframe")}
    if (name_column_groups_loading == PC_to_plot[1]) {stop(paste0('This is a problem..: name_column_groups_loading should not be named as "', PC_to_plot[1], '"'))}
    if (name_column_groups_loading == PC_to_plot[2]) {stop(paste0('This is a problem..: name_column_groups_loading should not be named as "', PC_to_plot[2], '"'))}

    if (!is.factor(pull(GCPlist$proteinINFO, name_column_groups_loading))) {
      GCPlist$proteinINFO[,name_column_groups_loading] <- as.factor(pull(GCPlist$proteinINFO, name_column_groups_loading))
    }

    if (!is.null(col_pal_loading)) {
      if (!is.character(col_pal_loading)) {stop("col_pal_loading must be a character vector")}
      if (any(is.na(col_pal_loading))) {stop("col_pal_loading must not contain NAs")}

    } else {
      if (used_the_long_vector_of_colors) {
        col_pal_loading <- build_long_vector_of_colors()[(length(col_pal)+1):length(build_long_vector_of_colors())]
      } else {
        col_pal_loading <- build_long_vector_of_colors()
      }
    }

    if (length(col_pal_loading)<length(levels(pull(GCPlist$proteinINFO, name_column_groups_loading)))) {
      stop(paste0("There are ", length(levels(pull(GCPlist$proteinINFO, name_column_groups_loading))), " groups, and only ", length(col_pal_loading), " colors have been specified in col_pal_loading"))
    } else {
      col_pal_loading <- col_pal_loading[1:length(levels(pull(GCPlist$proteinINFO, name_column_groups_loading)))]
    }

    if (!all(are.colors(col_pal_loading))) {stop(paste0('col_pal_loading must contain valid colors. In particular, these are not: "', paste(col_pal_fv[!are.colors(col_pal_fv)], collapse = '", "')), '"')}
  }

  if (!is.null(name_column_labels)) {
    if (length(name_column_labels)!=1) {stop("name_column_labels must be NULL or a character of length 1")}
    if (!is.character(name_column_labels)) {stop("name_column_labels must be NULL or a character of length 1")}
    if (is.na(name_column_labels)) {stop("name_column_labels must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_labels)) != 1) {stop("The name passed in name_column_labels must be a name of a column of the sampleINFO dataframe")}
    if (name_column_labels == PC_to_plot[1]) {stop(paste0('This is a problem..: name_column_labels should not be named as "', PC_to_plot[1], '"'))}
    if (name_column_labels == PC_to_plot[2]) {stop(paste0('This is a problem..: name_column_labels should not be named as "', PC_to_plot[2], '"'))}
    if (name_column_labels == "allwiththis_label") {stop("Please, just don't pass 'allwiththis_label' to name_column_labels, thanks!")}
    if (name_column_labels == "thesearethesamplenamesused") {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_labels, thanks!")}
  }


  if (!is.null(name_column_labels_loading)) {
    if (length(name_column_labels_loading)!=1) {stop("name_column_labels_loading must be NULL or a character of length 1")}
    if (!is.character(name_column_labels_loading)) {stop("name_column_labels_loading must be NULL or a character of length 1")}
    if (is.na(name_column_labels_loading)) {stop("name_column_labels_loading must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$proteinINFO) == name_column_labels_loading)) != 1) {stop("The name passed in name_column_labels_loading must be a name of a column of the proteinINFO dataframe")}
    if (name_column_labels_loading == PC_to_plot[1]) {stop(paste0('This is a problem..: name_column_labels_loading should not be named as "', PC_to_plot[1], '"'))}
    if (name_column_labels_loading == PC_to_plot[2]) {stop(paste0('This is a problem..: name_column_labels_loading should not be named as "', PC_to_plot[2], '"'))}
  }

  if (length(ellipses_on_loading)!=1) {stop("ellipses_on_loading must be exclusively TRUE or FALSE")}
  if (!is.logical(ellipses_on_loading)) {stop("ellipses_on_loading must be exclusively TRUE or FALSE")}
  if (is.na(ellipses_on_loading)) {stop("ellipses_on_loading must be exclusively TRUE or FALSE")}



  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw, name_first_column = "thesearethesamplenamesused")
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ, name_first_column = "thesearethesamplenamesused")
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }

  df_intensities_wg <- df_intensities

  if (!is.null(name_column_groups)) {
    df_intensities_wg <- add_column(df_intensities_wg,
                                    allwiththis = factor(NA, levels = levels(pull(GCPlist$sampleINFO, name_column_groups))),
                                    .after = 1)
    colnames(df_intensities_wg)[2] <- name_column_groups

    if (colnames(df_intensities_wg)[1] == colnames(df_intensities_wg)[2]) {
      the_new_name <- paste0(colnames(df_intensities_wg)[2], "_bis")
      colnames(df_intensities_wg)[2] <- the_new_name
      name_column_groups <- the_new_name
    }

    for (i in 1:length(pull(df_intensities_wg, 1))) {
      df_intensities_wg[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
    }
  }

  if (!is.null(name_column_labels)) {
    if (all(pull(df_intensities_wg, 1) == pull(GCPlist$sampleINFO, name_column_labels))) {
      colnames(df_intensities_wg)[1] <- name_column_labels

      if (colnames(df_intensities_wg)[1] == colnames(df_intensities_wg)[2]) {
        thiis_new_name <- paste0(colnames(df_intensities_wg)[1], "_bis")
        colnames(df_intensities_wg)[1] <- thiis_new_name
        name_column_labels <- thiis_new_name
      }

    } else {
      df_intensities_wg <- add_column(df_intensities_wg,
                                      allwiththis_label = NA,
                                      .after = 1)
      colnames(df_intensities_wg)[2] <- name_column_labels

      if (colnames(df_intensities_wg)[1] == colnames(df_intensities_wg)[2]) {
        this_new_name <- paste0(colnames(df_intensities_wg)[2], "_bis")
        colnames(df_intensities_wg)[2] <- this_new_name
        name_column_labels <- this_new_name
      }

      for (i in 1:length(pull(df_intensities_wg, 1))) {
        df_intensities_wg[i, name_column_labels] <- pull(GCPlist$sampleINFO, name_column_labels)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
      }
    }
  }


  PCA_list <- suppressWarnings(GetFeatistics::getPCA(df = df_intensities_wg,
                                                     v = colnames(df_intensities)[-1],
                                                     s = name_column_labels,
                                                     f = name_column_groups,
                                                     dfv = GCPlist$proteinINFO,
                                                     sv = name_column_labels_loading,
                                                     fv = name_column_groups_loading,
                                                     labels_on_loading = FALSE,
                                                     center = center,
                                                     scale. = scale.,
                                                     col_pal = col_pal,
                                                     col_pal_fv = col_pal_loading,
                                                     PC_to_plot = PC_to_plot,
                                                     ellipses_on_score = ellipses_on_score,
                                                     ellipses_on_loading = ellipses_on_loading))

  object_to_return <- NULL

  if (scores_or_loadings == "scores") {
    object_to_return <- PCA_list$score_plot
  } else if (scores_or_loadings == "loadings") {
    object_to_return <- PCA_list$loading_plot
  } else {
    stop('scores_or_loadings must be "scores" or "loadings"')
  }


  for (nlay in 1:length(object_to_return$layers)) {
    if (!inherits(object_to_return$layers[[nlay]]$geom, "GeomPoint")) {
      object_to_return$layers[[nlay]]$show.legend <- FALSE
    }
  }

  if (scores_or_loadings == "scores") {
    is_point <- vapply(object_to_return$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))
    object_to_return$layers[[which(is_point)[1]]]$aes_params$size <- 4
  }

  is_repel <- rep(FALSE, length(object_to_return$layers))
  for (i in seq_along(object_to_return$layers)) {
    is_repel[i] <- inherits(object_to_return$layers[[i]]$geom, "GeomTextRepel")
  }
  idx_repel <- which(is_repel)

  if (length(idx_repel)>0) {
    for (idx in idx_repel) {
      ly <- object_to_return$layers[[idx]]

      if (!is.null(ly$mapping) && "colour" %in% names(ly$mapping)) ly$mapping$colour <- NULL
      if (!is.null(ly$mapping) && "color"  %in% names(ly$mapping)) ly$mapping$color  <- NULL

      if (is.null(ly$mapping)) {ly$mapping <- ggplot2::aes()}
      ly$mapping$colour <- I("black")

      ly$aes_params$colour <- "black"

      object_to_return$layers[[idx]] <- ly
    }
  }

  object_to_return <- object_to_return +
    theme(panel.border = element_rect(linewidth = 1.2))


  return(object_to_return)
}
