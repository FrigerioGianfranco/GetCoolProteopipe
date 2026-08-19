# Filtering proteins

It filters the intensities data frame considering certain criteria.

## Usage

``` r
GCP_FilterProteins(GCPlist, operation)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- operation:

  character. Operation(s) to be applied considering a column of
  proteinINFO (for example: "ttest_PvaluesFDR \< 0.05"). All the TRUE
  from the operation will be kept.

## Value

The GCPlist with a potentially reduced number of rows in the intensities
tables.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist14f <- GCP_FilterProteins(GCPlist = GCPlist14,
                                 operation = "ttest_Pvalues < 0.05")


GCPlist14f2 <- GCP_FilterProteins(GCPlist = GCPlist14,
                                  operation = "`Gene names` %in% c('Ighg1', 'Ighm', 'Cdv3', 'Tceb1', 'Ktn1')")


} # }


```
