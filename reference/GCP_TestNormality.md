# Testing the normality

It performs the Shapiro-Wilk test and print out a density plot and a
qqplot.

## Usage

``` r
GCP_TestNormality(
  GCPlist,
  print_DensityPlot = FALSE,
  print_QQPlot = FALSE,
  print_only_these_protid = NULL,
  print_only_the_first = 10
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- print_DensityPlot:

  logical. If TRUE, a PDF with the density plots will be created in the
  current working directory.

- print_QQPlot:

  logical. If TRUE, a PDF with the Q-Q plots will be created in the
  current working directory.

- print_only_these_protid:

  NULL or a character vector of protid. Only the density/q-q-plot of the
  specified protid will be exported in the pdf files.

- print_only_the_first:

  NULL or a numeric integer. Only the density/q-q-plot of the those
  first proteins will be exported in the pdf files (suggested as usually
  with more than a thousand proteins, the file will be too big).

## Value

The GCPlist with the results of the Shapiro-Wilk test in the proteinINFO
table.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist05 <- GCP_TestNormality(GCPlist04)

} # }

```
