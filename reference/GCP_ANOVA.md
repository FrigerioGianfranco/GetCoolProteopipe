# performing an ANOVA on data

It performs a one-way analysis of variance on protein intensities, with
also a post-hoc Tukey's Test.

## Usage

``` r
GCP_ANOVA(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  FDR = TRUE,
  pcutoff = 0.05
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- name_column_groups:

  character of length 1. The name of the column of the sampleINFO table
  containing the sample groups. Since this is an ANOVA, there must be at
  least three groups.

- FDR:

  logical. If TRUE, after performing the ANOVA, it also correct p-values
  across the different proteins with a false discovery rate multiple
  comparison correction (method "fdr" of the function p.adjust).

- pcutoff:

  a numeric of length 1, must be between 0 and 1. The difference between
  paired groups will be reported only if the p-values is below the
  cut-off reported here.

## Value

The GCPlist with the results of the ANOVA added to the proteinINFO data
frame.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist14a1 <- GCP_ANOVA(GCPlist = GCPlist14a,
                         name_column_groups = "Group_multi")

} # }


```
