# Produce bar plot graphs out of Gene Ontologies (GO) results.

It creates bar plot graphs considering the output of the GCP_GO()
function.

## Usage

``` r
GCP_GO_barplot(
  GO_list,
  barfill = "pvalue",
  minCount = 1,
  maxDescriptions = Inf
)
```

## Arguments

- GO_list:

  A list obtained with the function GCP_GO().

- barfill:

  A name of a numeric column of the GO result data frames. Ideally one
  of: "pvalue", "p.adjust", "qvalue"

- minCount:

  numeric. A minimum number of Count to be included in the graph

- maxDescriptions:

  Numeric. Maximum number of bars to plot. rows with the lowest values
  of the column indicated in barfill and, in case of tie, with the
  higher number of Count, will be prioritised.

## Value

A list of ggplot objects.

## Examples

``` r

if (FALSE) { # \dontrun{
this_GO_list <- GCP_GO(GCPlist = GCPlist12a2,
                       name_column_groups = "Group_multi",
                       OrgDb = "mouse")

GO_bar_plots <- GCP_GO_barplot(GO_list = this_GO_list,
                               barfill = "pvalue")

export_figures(GO_bar_plots)

} # }

```
