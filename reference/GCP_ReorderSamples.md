# Reorder samples

Reorder samples in the CGP list.

## Usage

``` r
GCP_ReorderSamples(
  GCPlist,
  sample_names_ordered = NULL,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups")
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- sample_names_ordered:

  NULL or a character or factor vector containing the name of samples
  with the new desired order.

- name_column_groups:

  NULL or a character of length 1. The name of the column of the
  sampleINFO table containing the sample groups. If sample_names_ordered
  is NULL, the order will be based on the on such groups.

## Value

the GCPlist with samples reordered.

## Examples

``` r
if (FALSE) { # \dontrun{

# Reordering samples based on the groups (by default based on the set name_column_groups):

GCPlist02rd1 <- GCP_ReorderSamples(GCPlist02)

# Reordering samples specifying the samples:

GCPlist02rd2 <- GCP_ReorderSamples(GCPlist02, sample_names_ordered = c("S_1", "V_1", "S_2", "V_2", "S_3", "V_3", "S_4", "V_4", "S_5", "V_5"))


} # }

```
