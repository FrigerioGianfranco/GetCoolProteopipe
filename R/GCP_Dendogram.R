#' Get a dendogram of sample groups.
#'
#' It create a dendogram from hierarchical clustering of the data.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups.
#' @param col_pal NULL or a character vector containing colors. If NULL, colors from the pals package will be used (see function build_long_vector_of_colors).
#' @param rotate_names logical. If  TRUE, the names will be rotated vertically.
#'
#' @return A ggplot object.
#'
#'
#' @examples
#' \dontrun{
#'
#' Fig04_Dendogram_before_processing <- GCP_Dendogram(GCPlist05)
#' export_figures(Fig04_Dendogram_before_processing)
#'
#' }
#'
#'
#'
#' @import ggdendro
#'
#' @export
GCP_Dendogram <- function(GCPlist, name_column_groups = getOption("GetCoolProteopipe.name_column_groups"), col_pal = getOption("GetCoolProteopipe.col_pal"), rotate_names = TRUE) {

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
    if (name_column_groups == "allwiththis") {
      col_pal <- "black"
    } else {
      col_pal <- build_long_vector_of_colors()
    }
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

  if (length(rotate_names)!=1) {stop("rotate_names must be exclusively TRUE or FALSE")}
  if (!is.logical(rotate_names)) {stop("rotate_names must be exclusively TRUE or FALSE")}
  if (is.na(rotate_names)) {stop("rotate_names must be exclusively TRUE or FALSE")}

  df_intensities <- GetFeatistics::transpose_feat_table(GCPlist$intensities, name_first_column = colnames(GCPlist$sampleINFO)[1])


  df_intensities_wrn <- as.data.frame(df_intensities)
  rownames(df_intensities_wrn) <- pull(df_intensities, colnames(GCPlist$sampleINFO)[1])
  df_intensities_wrn <- df_intensities_wrn[, which(colnames(df_intensities_wrn) != colnames(GCPlist$sampleINFO)[1])]


  the_model <- hclust(dist(df_intensities_wrn), "ave")
  the_dhc <- as.dendrogram(the_model)
  the_dendro_data <- dendro_data(the_dhc, type = "rectangle")
  the_dendro_data$labels[, name_column_groups] <- factor(rep(NA, length(pull(the_dendro_data$labels, 1))), levels = levels(pull(GCPlist$sampleINFO, name_column_groups)))

  for (i in 1:length(pull(the_dendro_data$labels, 1))) {
    the_dendro_data$labels[i, name_column_groups] <- pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1) == the_dendro_data$labels$label[i])]
  }

  the_plot <- ggplot(segment(the_dendro_data)) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
    theme_dendro()

  if (rotate_names) {
    the_plot <- the_plot +
      geom_text(data = the_dendro_data$labels,
                aes(x = x, y = y, label = label, color = !!sym(name_column_groups)),
                hjust = 0, vjust = 0, angle = 270, show.legend = FALSE) +
      scale_color_manual(values = named_colors_bygroup) +
      scale_y_continuous(expand = expansion(mult = c(0.014*max(nchar(the_dendro_data$labels$label)), 0)))

  } else {
    the_plot <- the_plot +
      geom_text(data = the_dendro_data$labels,
                aes(x = x, y = y, label = label, colour = !!sym(name_column_groups)),
                hjust = 0.5, vjust = 1, angle = 0, show.legend = FALSE) +
      scale_color_manual(values = named_colors_bygroup)
  }

  if (name_column_groups != "allwiththis") {

    the_plot <- the_plot +
      geom_point(data = the_dendro_data$labels, aes(x = x, y = y, color = !!sym(name_column_groups)), alpha = 0, show.legend = TRUE) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, size = 3)))
  }

  return(the_plot)
}


