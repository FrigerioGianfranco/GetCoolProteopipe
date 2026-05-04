# Get ggven for protein overlapping among groups.

It creates venn diagrams to show the overlapping proteins per sample
group.

## Usage

``` r
GCP_Venn(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  col_pal = getOption("GetCoolProteopipe.col_pal"),
  consider_genenames = FALSE,
  auto_scale_circles = FALSE,
  add_title = FALSE
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- name_column_groups:

  character of length 1. The name of the column of the sampleINFO table
  containing the sample groups. The sample groups must be between 2 and
  4.

- col_pal:

  NULL or a character vector containing colors. If NULL, colors from the
  pals package will be used (see function build_long_vector_of_colors).

- consider_genenames:

  logical. If TRUE, the gene names will be considered to perform the
  grouping, instead of the protid.

- auto_scale_circles:

  logical. If TRUE and if there are only 2 groups, the dimensions of
  circles will be scaled to the number of intersections.

- add_title:

  logical. If TRUE, a short default title is added at the top of the
  graph.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{

Fig08_Venn_post_filteringNA <- GCP_Venn(GCPlist07)
export_figures(Fig08_Venn_post_filteringNA)

} # }


```
