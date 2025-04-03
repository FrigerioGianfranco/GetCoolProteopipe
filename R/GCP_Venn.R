#' Get ggven for protein overlapping among groups.
#'
#' It creates venn diagrams to show the overlapping proteins per sample group.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups character of length 1. The name of the column of the sampleINFO table containing the sample groups. The sample groups must be between 2 and 4.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The venn diagram will be performed only in the specified data.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used (see function build_long_vector_of_colors).
#' @param auto_scale_circles logical. If TRUE and if there are only 2 groups, the dimensions of circles will be scaled to the number of intersections.
#'
#' @return A ggplot object.
#'
#' @import ggvenn
#'
#' @export
GCP_Venn <- function(GCPlist, name_column_groups = NULL, raw_or_LFQ = c("raw", "LFQ"), col_pal = NULL, auto_scale_circles = FALSE) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
    if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
    if (is.na(name_column_groups)) {stop("name_column_groups must be a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

    if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
      GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
    }

    if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) < 2 | length(levels(pull(GCPlist$sampleINFO, name_column_groups))) > 4) {stop("The sample groups must be only 2, 3 or 4")}

  } else {
    stop("please, specifiy the name of the sample group column in the name_column_groups argoument")
  }

  if (!identical(tolower(raw_or_LFQ), c("raw", "lfq"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
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

  if (length(auto_scale_circles)!=1) {stop("auto_scale_circles must be exclusively TRUE or FALSE")}
  if (!is.logical(auto_scale_circles)) {stop("auto_scale_circles must be exclusively TRUE or FALSE")}
  if (is.na(auto_scale_circles)) {stop("auto_scale_circles must be exclusively TRUE or FALSE")}


  if (raw_or_LFQ == "raw") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_raw, name_first_column = colnames(GCPlist$sampleINFO)[1])
  } else if (raw_or_LFQ == "lfq") {
    df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$quant_LFQ, name_first_column = colnames(GCPlist$sampleINFO)[1])
  } else {
    stop('raw_or_LFQ must be one of "raw", "LFQ"')
  }

  df_intensities_w_groups <- left_join(x = GCPlist$sampleINFO[,c(colnames(GCPlist$sampleINFO)[1], name_column_groups)], y = df_intensities, by = colnames(GCPlist$sampleINFO)[1], suffix = c("_INFO", "_intensities"))

  list_for_ggven <- vector(mode = "list", length = length(levels(pull(GCPlist$sampleINFO, name_column_groups))))
  names(list_for_ggven) <- levels(pull(GCPlist$sampleINFO, name_column_groups))

  for (gr in levels(pull(GCPlist$sampleINFO, name_column_groups))) {

    df_intensities_w_groups_fill <- df_intensities_w_groups[which(pull(df_intensities_w_groups, name_column_groups) == gr), colnames(df_intensities_w_groups)[which(!colnames(df_intensities_w_groups)%in%c(colnames(GCPlist$sampleINFO)[1], name_column_groups))]]

    list_for_ggven[[gr]] <- character()

    for (a in colnames(df_intensities_w_groups_fill)) {
      if (any(!is.na(pull(df_intensities_w_groups_fill, a)))) {
        list_for_ggven[[gr]] <- c(list_for_ggven[[gr]], a)
      }
    }
  }


  fill_vector <- col_pal
  names(fill_vector) <- NULL




  if (length(levels(pull(GCPlist$sampleINFO, name_column_groups))) == 2 & auto_scale_circles) {
    ggven_graph <- ggvenn(list_for_ggven, fill_color = fill_vector, auto_scale = TRUE) + theme_bw() + theme(axis.title = element_blank(),
                                                                                                            axis.text = element_blank(),
                                                                                                            axis.ticks = element_blank(),
                                                                                                            axis.line = element_blank(),
                                                                                                            panel.grid.major = element_blank(),
                                                                                                            panel.grid.minor = element_blank(),
                                                                                                            panel.border = element_blank(),
                                                                                                            plot.background = element_blank())
  } else {
    ggven_graph <- ggvenn(list_for_ggven, fill_color = fill_vector) + theme_bw() + theme(axis.title = element_blank(),
                                                                                         axis.text = element_blank(),
                                                                                         axis.ticks = element_blank(),
                                                                                         axis.line = element_blank(),
                                                                                         panel.grid.major = element_blank(),
                                                                                         panel.grid.minor = element_blank(),
                                                                                         panel.border = element_blank(),
                                                                                         plot.background = element_blank())
  }


  if (raw_or_LFQ == "raw") {

    ggven_graph <- ggven_graph + ggtitle("Proteins in groups, raw data") + theme(plot.title = element_text(hjust = 0.5))

  } else if (raw_or_LFQ == "lfq") {

    ggven_graph <- ggven_graph + ggtitle("Proteins in groups, LFQ data") + theme(plot.title = element_text(hjust = 0.5))
  }


  return(ggven_graph)
}
