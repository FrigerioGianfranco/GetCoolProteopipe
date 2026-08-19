# performing a t-test on data

It performs a t-test analyses on the proteins intensities.

## Usage

``` r
GCP_ttest(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  paired = FALSE,
  FDR = TRUE,
  pcutoff = 0.05
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- name_column_groups:

  character of length 1. The name of the column of the sampleINFO table
  containing the sample groups. Since this is a t-test, there must be
  exactly two groups.

- paired:

  logical. If FALSE it performs non-paired t-tests. If TRUE it performs
  paired t-tests.

- FDR:

  logical. If TRUE, after performing the t-tests, it also correct
  p-values across the different proteins with a false discovery rate
  multiple comparison correction (method "fdr" of the function
  p.adjust).

- pcutoff:

  a numeric of length 1, must be between 0 and 1. The difference between
  groups will be reported only if the p-values is below the cut-off
  reported here.

## Value

The GCPlist with the results of the t-test added to the proteinINFO data
frame.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist13 <- GCP_ttest(GCPlist = GCPlist12,
                       paired = TRUE,
                       FDR = TRUE,
                       pcutoff = 0.05)

} # }


```
