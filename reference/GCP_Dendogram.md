# Get a dendogram of sample groups.

It create a dendogram from hierarchical clustering of the data.

## Usage

``` r
GCP_Dendogram(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  col_pal = getOption("GetCoolProteopipe.col_pal"),
  rotate_names = TRUE
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups.

- col_pal:

  NULL or a character vector containing colors. If NULL, colors from the
  pals package will be used (see function build_long_vector_of_colors).

- rotate_names:

  logical. If TRUE, the names will be rotated vertically.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{

Fig04_Dendogram_before_processing <- GCP_Dendogram(GCPlist05)
export_figures(Fig04_Dendogram_before_processing)

} # }


```
