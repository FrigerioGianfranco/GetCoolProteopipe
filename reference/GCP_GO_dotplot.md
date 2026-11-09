# Produce dot plot graphs out of Gene Ontologies (GO) results.

It creates dot plot graphs considering the output of the GCP_GO()
function.

## Usage

``` r
GCP_GO_dotplot(
  GO_list,
  dotcol = "pvalue",
  xaxis = "GeneRatio",
  differential_together = TRUE,
  minCount = 1,
  maxDescriptions = Inf
)
```

## Arguments

- GO_list:

  A list obtained with the function GCP_GO().

- dotcol:

  A name of a numeric column of the GO result data frames. Ideally one
  of: "pvalue", "p.adjust", "qvalue".

- xaxis:

  A name of a column of the GO result data frames, such as "GeneRatio"
  or "BgRatio".

- differential_together:

  logical. If TRUE, and if the results allow it, the differential
  analyses of the same comparison will be considered togheter in the
  same graph.

- minCount:

  numeric. A minimum number of Count to be included in the graph.

- maxDescriptions:

  Numeric. Maximum number of lines to plot. rows with the lowest values
  of the column indicated in dotcol and, in case of tie, with the higher
  number of Count, will be prioritised.

## Value

A list of ggplot objects.

## Examples

``` r

if (FALSE) { # \dontrun{
this_GO_list <- GCP_GO(GCPlist = GCPlist12a2,
                       name_column_groups = "Group_multi",
                       OrgDb = "mouse")

GO_dot_plots <- GCP_GO_dotplot(GO_list = this_GO_list,
                               dotcol = "pvalue",
                               xaxis = "GeneRatio",
                               differential_together = TRUE)

export_figures(GO_dot_plots)

} # }

```
