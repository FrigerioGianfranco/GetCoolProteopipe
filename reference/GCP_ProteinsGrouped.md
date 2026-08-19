# Know proteins shared by groups.

It creates an additional column in the proteinINFO data frame containing
where in which groups proteins are present.

## Usage

``` r
GCP_ProteinsGrouped(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups")
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- name_column_groups:

  character of length 1. The name of the column of the sampleINFO table
  containing the sample groups. The sample groups must be between at
  least 2.

## Value

The GCPlist with an additional column in the proteinINFO data frame.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist08 <- GCP_ProteinsGrouped(GCPlist07)

} # }


```
