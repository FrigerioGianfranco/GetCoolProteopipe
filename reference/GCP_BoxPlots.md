# Get Box plots from intensities.

It creates box plots showing the distribution of proteins in each group.

## Usage

``` r
GCP_BoxPlots(
  GCPlist,
  by_samples = TRUE,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  col_pal = getOption("GetCoolProteopipe.col_pal"),
  Title = "Distribution of intensities",
  only_these_protid = NULL,
  only_the_first = NULL
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- by_samples:

  logical. If TRUE, a single box plot graph will be returned with a box
  for each sample. If FALSE, multiple box plots can be generated for
  each protein, as specified in the only_these_protid or only_the_first
  arguments.

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups.

- col_pal:

  NULL or a character vector containing colors. If NULL, colors from the
  pals package will be used. To see the colors, run
  build_long_vector_of_colors().

- Title:

  NULL or character of length 1. The title you want to ad on the top of
  the plot. Used only if by_samples is TRUE.

- only_these_protid:

  NULL or a character vector of protid. Box plots will be performed only
  for those proteins. Used only if by_samples is FALSE.

- only_the_first:

  NULL or a numeric integer. If only_these_protid is NULL, box plots
  will be performed only for those first proteins. Used only if
  by_samples is FALSE.

## Value

A ggplot object, if by_samples is TRUE; or a list of ggplot objects, if
by_samples is FALSE.

## Examples

``` r
if (FALSE) { # \dontrun{

# A single Box-plot graph with distribution of all proteins for each sample:

Fig03_BoxPlot_before_processing <- GCP_BoxPlots(GCPlist = GCPlist05,
                                                by_samples = TRUE,
                                                Title = "Distribution of intensities, before processing")
export_figures(Fig03_BoxPlot_before_processing)


# Multiple Box-plot graphs, each with the distributions of a single protein with boxes by sample groups:

Fig17_BoxPlots_sign <- GCP_BoxPlots(GCPlist = GCPlist14f,
                                    by_samples = FALSE)
export_figures(Fig17_BoxPlots_sign)                            ## exporting each box-plot as a single png file
export_figures(Fig17_BoxPlots_sign, exprt_fig_type = "pdf")    ## exporting all box-plots in a pdf file with many pages

} # }


```
