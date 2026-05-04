# Imputation of missing values

It performs the imputation of missing values, using the same method of
the functions scImpute and tImpute from the PhosR package.

## Usage

``` r
GCP_NAimputation(
  GCPlist,
  quant_rate = 0.5,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  seed = 123,
  method = c("both", "scImpute", "tImpute")
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- quant_rate:

  numeric, between 0 an 1 (included). Quantification rate per group,
  considered for the scImpute function.

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups.

- seed:

  numeric. The value that will be used for set.seed.

- method:

  one of: "both", "scImpute", "tImpute". If both is selected, it
  performs first the scImpute, and then the tImpute.

## Value

a GCPlist list with the missing values imputed.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist09 <- GCP_NAimputation(GCPlist = GCPlist08,
                              quant_rate = 0.5)

} # }


```
