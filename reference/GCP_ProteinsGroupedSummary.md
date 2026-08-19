# Get a summary out after the GCP_ProteinsGrouped function.

It creates a small table with the summary of the results of the
GCP_ProteinsGrouped function. PLEASE NOTE that you must run the function
GCP_ProteinsGrouped before!

## Usage

``` r
GCP_ProteinsGroupedSummary(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  digits_perc = NULL,
  add_perc_symbol = FALSE
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- name_column_groups:

  character of length 1. The name of the column of the sampleINFO table
  containing the sample groups.

- digits_perc:

  NULL or numeric integer of length 1. The number of digits to round the
  percentages.

- add_perc_symbol:

  logical. If TRUE, it will add the ' %' symbol at the percentages.

## Value

A tibble with summary of the presence of proteins in indicated groups.

## Examples

``` r
if (FALSE) { # \dontrun{

SummaryTable_ProteinsGrouped_filteringNA <- GCP_ProteinsGroupedSummary(GCPlist08)
export_the_table(SummaryTable_ProteinsGrouped_filteringNA)

} # }


```
