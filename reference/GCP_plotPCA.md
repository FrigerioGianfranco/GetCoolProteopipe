# Plot Principal Component analysis.

It plots a principal component analysis.

## Usage

``` r
GCP_plotPCA(
  GCPlist,
  scores_or_loadings = c("scores", "loadings"),
  PC_to_plot = c("PC1", "PC2"),
  center = TRUE,
  scale. = FALSE,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  name_column_labels = NULL,
  col_pal = getOption("GetCoolProteopipe.col_pal"),
  ellipses_on_score = TRUE,
  name_column_groups_loading = NULL,
  name_column_labels_loading = NULL,
  col_pal_loading = NULL,
  ellipses_on_loading = FALSE
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- scores_or_loadings:

  one of the following: "scores", "loadings". Specify here if you want
  to plot the scores or the loadings.

- PC_to_plot:

  character of length 2. Specify here the two principal components to
  plot.

- center:

  logical. Whether the variables should be shifted to be zero centered
  (as in the prcomp function).

- scale.:

  logical. whether the variables should be scaled to have unit variance
  before the analysis takes place (as in prcomp function).

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups. Specify it only if you
  want to color the points of the score plot.

- name_column_labels:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample names. Specify it only if you
  want to add a label to the points of the score plot.

- col_pal:

  a character vector containing colors for groups for the score plot. If
  NULL, colors from the pals package will be used (see function
  GetFeatistics::build_long_vector_of_colors).

- ellipses_on_score:

  logical. If you specified name_column_groups and this is TRUE,
  ellipses will be added to the score plot.

- name_column_groups_loading:

  NULL or character of length 1. The name of the column of the
  proteinINFO table containing the protein groups. Specify it only if
  you want to color the points of the loading plot.

- name_column_labels_loading:

  NULL or character of length 1. The name of the column of the
  proteinINFO table containing the protein/protein names. Specify it
  only if you want to add a label to the points of the loading plot.

- col_pal_loading:

  a character vector containing colors for groups for the loading plot.
  If NULL, colors from the pals package will be used (see function
  GetFeatistics::build_long_vector_of_colors).

- ellipses_on_loading:

  logical. If you specified name_column_groups_loading and this is TRUE,
  ellipses will be added to the loading plot.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{

# to get a score plot:

Fig13_PCA_score_plot <- GCP_plotPCA(GCPlist = GCPlist11,
                                    scores_or_loadings = "scores",
                                    PC_to_plot = c("PC1", "PC2"),
                                    name_column_labels = "Sample",
                                    ellipses_on_score = TRUE)
export_figures(Fig13_PCA_score_plot)


# to get a loading plot:

Fig14_PCA_loading_plot <- GCP_plotPCA(GCPlist = GCPlist11,
                                      scores_or_loadings = "loadings",
                                      PC_to_plot = c("PC1", "PC2"),
                                      name_column_groups_loading = NULL,
                                      name_column_labels_loading = "Protein names",
                                      col_pal_loading = NULL,
                                      ellipses_on_loading = FALSE)
export_figures(Fig14_PCA_loading_plot)

} # }


```
