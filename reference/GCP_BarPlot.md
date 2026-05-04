# Get bar plots with protein numbers.

It create bar plots to show the number of valid proteins per samples.

## Usage

``` r
GCP_BarPlot(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  col_pal = getOption("GetCoolProteopipe.col_pal"),
  bar_width = NULL,
  label_numbers = TRUE,
  label_numbers_size = NULL,
  showCV = FALSE,
  rotate_sample_names = FALSE
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups.

- col_pal:

  NULL or a character vector containing colors. If NULL, colors from the
  pals package will be used (see function build_long_vector_of_colors).

- bar_width:

  NULL or a number. You can pass here the width of the bars.

- label_numbers:

  logical. If TRUE, it adds the number of proteins on the top of each
  bar.

- label_numbers_size:

  NULL or numeric of length 1. If specified and if label_numbers is
  TRUE, this is the size of the numbers of proteins of the top of each
  bar.

- showCV:

  logical. If TRUE, it adds to the graphs the coefficient of variation
  (CV%) of the number of proteins for each group.

- rotate_sample_names:

  logical. If TRUE, it rotate sample name labels of x-axis in vertical
  position.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{

Fig07_BarPlot_post_filteringNA <- GCP_BarPlot(GCPlist = GCPlist07,
                                              label_numbers = TRUE,
                                              showCV = TRUE,
                                              rotate_sample_names = TRUE)
export_figures(Fig07_BarPlot_post_filteringNA)

} # }


```
