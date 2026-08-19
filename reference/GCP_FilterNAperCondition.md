# Filter proteins based on missing values per condition

In the intensities table, it keeps only proteins that are not missing
values for a defined ratio in at least one of the specified sample
groups.

## Usage

``` r
GCP_FilterNAperCondition(
  GCPlist,
  ratio = 0.5,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups")
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- ratio:

  numeric, between 0 an 1 (included). Ratio of non-missing values, per
  sample group, wanted to pass this filtration step (the higher this
  ratio, the less proteins will be kept).

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups.

## Value

a GCPlist list with filtered tables.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist07 <- GCP_FilterNAperCondition(GCPlist = GCPlist06,
                                      ratio = 0.5)

} # }


```
