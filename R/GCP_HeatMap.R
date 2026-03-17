#' Compute a Heat Map.
#'
#' It performs a heat map out of GCPdata.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character vector. Name(s) of the column(s) of the sampleINFO table. Specify it if you want to have those groups specified on the heatmap.
#' @param name_column_labels NULL or a character of length 1. The name of the column of sampleINFO containing the sample names. You need to pass it only if you want sample names on the heat map.
#' @param name_column_groups_protein NULL or character vector. Name(s) of the column(s) of the proteinINFO table. Specify it if you want to have those groups specified on the heatmap.
#' @param name_column_labels_protein NULL or a character of length 1. The name of the column of proteinINFO containing the protein names. You need to pass it only if you want sample names on the heat map.
#' @param order_samples_by NULL or character. The name(s) of the column(s) of sampleINFO that you want to use to order the values.
#' @param order_protein_by NULL or character. The name(s) of the column(s) of proteinINFO that you want to use to order the values.
#' @param trnsp logical. If TRUE, proteins will be in rows and samples in columns of the heat map; if FALSE, the opposite.
#' @param cluster_rows logical. If TRUE, a dendrogram built from the hierarchical clustering of distances of row values will be added to the left side of the heat map.
#' @param cluster_columns logical. If TRUE, a dendrogram built from the hierarchical clustering of distances of column values will be added above the heat map.
#' @param name_rows logical. If TRUE, names will be added to rows, on the right of the heat map.
#' @param name_columns logical. if TRUE, names will be added to column, below the heat map.
#' @param rotate_name_columns  logical. if TRUE, column names will be rotated vertically (this argument is meaningless if name_columns is FALSE).
#' @param three_heat_colors character of length 3, each specifying a color. These 3 colors will be used as color scale for values of the heat map.
#' @param set_heat_colors_limits logical. If TRUE, the absolute of the minimum or the absolute of the maximum value (which is higher) will be set in positive as the upper limit and in negative as the lower limit for the color gradients of values of the heat map (this will also set the middle color exactly to zero).
#' @param heat_colors_limits NULL or a numeric of length 2. If set_heat_colors_limits is FALSE, you can specify here the limits for the color gradients of values of the heat map (if NULL, the maximum and the minimum values will be used).
#' @param col_pal_list NULL, a character vector, or a list of character vector. If a list, elements should be named as name_column_groups and name_column_groups_protein; each element has to be a character vector containing colors. Those colors will be used for the rectangles of the group classifications. For each element, if NULL, colors will be taken from the pals package (see the function build_long_vector_of_colors).
#'
#' @return a ggplot object.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#' Fig16_The_Heat_Map_sign <- GCP_HeatMap(GCPlist = GCPlist14f,
#'                                        name_column_groups = c("Condition", "Replicate"),
#'                                        name_column_labels = "Sample",
#'                                        name_column_groups_protein = NULL,
#'                                        name_column_labels_protein = "Protein names",
#'                                        name_rows = FALSE,
#'                                        name_columns = TRUE,
#'                                        rotate_name_columns = TRUE,
#'                                        col_pal_list = list(Condition = c(S = "green", V = "blue"),
#'                                                            Replicate = c("lightskyblue", "lightskyblue1", "lightskyblue2", "lightskyblue3", "lightskyblue4")))
#' export_figures(Fig16_The_Heat_Map_sign)
#'
#' }
#'
#'
#'
#' @export
GCP_HeatMap <- function(GCPlist,
                        name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), name_column_labels = NULL,
                        name_column_groups_protein = NULL, name_column_labels_protein = NULL,
                        order_samples_by = NULL, order_protein_by = NULL,
                        trnsp = TRUE, cluster_rows = TRUE, cluster_columns = TRUE,
                        name_rows = FALSE, name_columns = FALSE, rotate_name_columns = TRUE,
                        three_heat_colors = c("red", "white", "blue"), set_heat_colors_limits = FALSE, heat_colors_limits = NULL,
                        col_pal_list = getOption("GetCoolProteopipe.col_pal")) {

  checkGCPlist(GCPlist)

  if (!is.null(name_column_groups)) {
    if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character")}
    if (any(is.na(name_column_groups))) {stop("name_column_groups must not contain NAs")}
    if (length(name_column_groups)==1) {
      cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
    } else {
      cat(paste0("\n -- The name_column_groups considered is c('", paste0(name_column_groups, collapse = "', '"), "') --\n\n"))
    }
    if (!all(name_column_groups %in% colnames(GCPlist$sampleINFO))) {stop("The names passed in name_column_groups must names of columns of the sampleINFO dataframe")}
    if (any(duplicated(colnames(GCPlist$sampleINFO)[which(colnames(GCPlist$sampleINFO)%in%name_column_groups)]))) {stop("there are duplicates in the name of sampleINFO considering the name_column_groups")}
    if (any(name_column_groups == "allwiththis")) {stop("Please, just don't pass 'allwiththis' to name_column_groups, thanks!")}
    if (any(name_column_groups == "thesearethesamplenamesused")) {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_groups, thanks!")}
  } else {
    cat("\n -- The name_column_groups considered is NULL --\n\n")
  }

  if (!is.null(name_column_labels)) {
    if (length(name_column_labels)!=1) {stop("name_column_labels must be NULL or a character of length 1")}
    if (!is.character(name_column_labels)) {stop("name_column_labels must be NULL or a character of length 1")}
    if (is.na(name_column_labels)) {stop("name_column_labels must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$sampleINFO) == name_column_labels)) != 1) {stop("The name passed in name_column_labels must be a name of a column of the sampleINFO dataframe")}
    if (name_column_labels == "allwiththis_label") {stop("Please, just don't pass 'allwiththis_label' to name_column_labels, thanks!")}
    if (name_column_labels == "thesearethesamplenamesused") {stop("Please, just don't pass 'thesearethesamplenamesused' to name_column_labels, thanks!")}
  }

  if (!is.null(name_column_groups_protein)) {
    if (!is.character(name_column_groups_protein)) {stop("name_column_groups_protein must be NULL or a character")}
    if (any(is.na(name_column_groups_protein))) {stop("name_column_groups_protein must not contain NAs")}
    if (!all(name_column_groups_protein %in% colnames(GCPlist$proteinINFO))) {stop("The names passed in name_column_groups_protein must names of columns of the proteinINFO dataframe")}
    if (any(duplicated(colnames(GCPlist$proteinINFO)[which(colnames(GCPlist$proteinINFO)%in%name_column_groups_protein)]))) {stop("there are duplicates in the name of proteinINFO considering the name_column_groups_protein")}
  }

  if (!is.null(name_column_labels_protein)) {
    if (length(name_column_labels_protein)!=1) {stop("name_column_labels_protein must be NULL or a character of length 1")}
    if (!is.character(name_column_labels_protein)) {stop("name_column_labels_protein must be NULL or a character of length 1")}
    if (is.na(name_column_labels_protein)) {stop("name_column_labels_protein must be NULL or a character of length 1, not a NA")}
    if (length(which(colnames(GCPlist$proteinINFO) == name_column_labels_protein)) != 1) {stop("The name passed in name_column_labels_protein must be a name of a column of the sampleINFO dataframe")}
  }

  if (!is.null(order_samples_by)) {
    if (!is.character(order_samples_by)) {stop("if not NULL, order_samples_by must be a character")}
    if (any(is.na(order_samples_by))) {stop("if not NULL, order_samples_by must not contain mising values")}
    if (!all(order_samples_by %in% colnames(GCPlist$sampleINFO))) {stop("if not NULL, the names you indicate in order_samples_by must correspond to names of columns in sampleINFO")}
  }

  if (!is.null(order_protein_by)) {
    if (!is.character(order_protein_by)) {stop("if not NULL, order_protein_by must be a character")}
    if (any(is.na(order_protein_by))) {stop("if not NULL, order_protein_by must not contain mising values")}
    if (!all(order_protein_by %in% colnames(GCPlist$proteinINFO))) {stop("if not NULL, the names you indicate in order_protein_by must correspond to names of columns in proteinINFO")}
  }


  if (length(trnsp)!=1) {stop("trnsp must be exclusively TRUE or FALSE")}
  if (!is.logical(trnsp)) {stop("trnsp must be exclusively TRUE or FALSE")}
  if (is.na(trnsp)) {stop("trnsp must be exclusively TRUE or FALSE")}

  if (length(cluster_rows)!=1) {stop("cluster_rows must be exclusively TRUE or FALSE")}
  if (!is.logical(cluster_rows)) {stop("cluster_rows must be exclusively TRUE or FALSE")}
  if (is.na(cluster_rows)) {stop("cluster_rows must be exclusively TRUE or FALSE")}

  if (length(cluster_columns)!=1) {stop("cluster_columns must be exclusively TRUE or FALSE")}
  if (!is.logical(cluster_columns)) {stop("cluster_columns must be exclusively TRUE or FALSE")}
  if (is.na(cluster_columns)) {stop("cluster_columns must be exclusively TRUE or FALSE")}

  if (length(name_rows)!=1) {stop("name_rows must be exclusively TRUE or FALSE")}
  if (!is.logical(name_rows)) {stop("name_rows must be exclusively TRUE or FALSE")}
  if (is.na(name_rows)) {stop("name_rows must be exclusively TRUE or FALSE")}

  if (length(name_columns)!=1) {stop("name_columns must be exclusively TRUE or FALSE")}
  if (!is.logical(name_columns)) {stop("name_columns must be exclusively TRUE or FALSE")}
  if (is.na(name_columns)) {stop("name_columns must be exclusively TRUE or FALSE")}

  if (name_columns) {
    if (length(rotate_name_columns)!=1) {stop("rotate_name_columns must be exclusively TRUE or FALSE")}
    if (!is.logical(rotate_name_columns)) {stop("rotate_name_columns must be exclusively TRUE or FALSE")}
    if (is.na(rotate_name_columns)) {stop("rotate_name_columns must be exclusively TRUE or FALSE")}
  }

  are.colors <- function (vect) {
    map_lgl(vect, ~tryCatch({
      is.matrix(col2rgb(.)) & ncol(col2rgb(.))>0
    }, error = function(e) {
      FALSE
    }))
  }

  if (!is.character(three_heat_colors)) {stop("three_heat_colors must be a character vector containing three colors")}
  if (length(three_heat_colors) != 3) {stop("three_heat_colors must be a character vector containing three colors")}
  if (any(is.na(three_heat_colors))) {stop("three_heat_colors must be a character vector containing three colors, with no missing values")}
  if (!all(are.colors(three_heat_colors))) {stop(paste0('three_heat_colors must be a character vector containing three colors. The following are not: "',
                                                        paste0(three_heat_colors[which(!are.colors(three_heat_colors))], collapse = '", "'), '"'))}

  if (length(set_heat_colors_limits)!=1) {stop("set_heat_colors_limits must be exclusively TRUE or FALSE")}
  if (!is.logical(set_heat_colors_limits)) {stop("set_heat_colors_limits must be exclusively TRUE or FALSE")}
  if (is.na(set_heat_colors_limits)) {stop("set_heat_colors_limits must be exclusively TRUE or FALSE")}

  if (set_heat_colors_limits==FALSE) {
    if (!is.null(heat_colors_limits)) {
      if (length(heat_colors_limits) != 2) {stop("if not NULL, heat_colors_limits must be a numeric of length 2")}
      if (!is.numeric(heat_colors_limits)) {stop("if not NULL, heat_colors_limits must be a numeric of length 2")}
    }
  }

  if (!is.null(name_column_groups) | !is.null(name_column_groups_protein)) {
    if (!is.null(col_pal_list)) {
      if (!is.list(col_pal_list)) {
        if ((length(name_column_groups) + length(name_column_groups_protein))==1) {
          if (!is.null(name_column_groups)) {
            col_pal_list <- list(aa = col_pal_list)
            names(col_pal_list) <- name_column_groups
          } else if (!is.null(name_column_groups_protein)) {
            col_pal_list <- list(aa = col_pal_list)
            names(col_pal_list) <- name_column_groups_protein
          }
        } else {
          stop("if not NULL, col_pal_list must be a list with character vector of colors to match name_column_groups and name_column_groups_protein")
        }
      }
      cat("\n -- col_pal_list is a list of colors:\n")
      print(col_pal_list)
      cat("--\n\n")
    } else {
      cat("\n -- col_pal_list is  NULL --\n\n")
    }
  }


  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities, name_first_column = "thesearethesamplenamesused")
  used_protid <- GCPlist$intensities$protid

  df_intensities_wg <- df_intensities

  if (!is.null(name_column_groups)) {
    for (ncg in rev(name_column_groups)) {

      if (!is.factor(pull(GCPlist$sampleINFO, ncg))) {

        GCPlist$sampleINFO <- mutate_at(GCPlist$sampleINFO, ncg, as.factor)
      }

      df_intensities_wg <- add_column(df_intensities_wg,
                                      allwiththis = factor(NA_character_, levels = levels(pull(GCPlist$sampleINFO, ncg))),
                                      .after = 1)
      colnames(df_intensities_wg)[2] <- ncg

      for (i in 1:length(pull(df_intensities_wg, 1))) {
        df_intensities_wg[i, ncg] <- pull(GCPlist$sampleINFO, ncg)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
      }
    }
  }

  if (!is.null(name_column_labels)) {
    if (all(pull(df_intensities_wg, 1) == pull(GCPlist$sampleINFO, name_column_labels))) {
      colnames(df_intensities_wg)[1] <- name_column_labels
    } else {
      df_intensities_wg <- add_column(df_intensities_wg,
                                      allwiththis_label = NA,
                                      .after = 1)
      colnames(df_intensities_wg)[2] <- name_column_labels

      for (i in 1:length(pull(df_intensities_wg, 1))) {
        df_intensities_wg[i, name_column_labels] <- pull(GCPlist$sampleINFO, name_column_labels)[which(pull(GCPlist$sampleINFO, 1) == pull(df_intensities_wg, 1)[i])]
      }
    }
  }

  the_heat_map <- GetFeatistics::getHeatMap(df = df_intensities_wg,
                                            v = colnames(df_intensities)[-1],
                                            s = name_column_labels,
                                            f = name_column_groups,
                                            dfv = GCPlist$proteinINFO[which(GCPlist$proteinINFO$protid %in% used_protid),],
                                            sv = name_column_labels_protein,
                                            fv = name_column_groups_protein,
                                            order_df_by = order_samples_by,
                                            order_dfv_by = order_protein_by,
                                            trnsp = trnsp,
                                            cluster_rows = cluster_rows,
                                            cluster_columns = cluster_columns,
                                            name_rows = name_rows,
                                            name_columns = name_columns,
                                            rotate_name_columns = rotate_name_columns,
                                            three_heat_colors = three_heat_colors,
                                            set_heat_colors_limits = set_heat_colors_limits,
                                            heat_colors_limits = heat_colors_limits,
                                            col_pal_list = col_pal_list)

  return(the_heat_map)
}
