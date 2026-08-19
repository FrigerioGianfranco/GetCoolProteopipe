# Change the name of samples.

Change the name of some samples.

## Usage

``` r
GCP_ChangeSampleNames(
  GCPlist,
  old_names = pull(GCPlist$sampleINFO, 1),
  new_names = paste0(old_names, "_updated"),
  old_new_names_table = NULL
)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- old_names:

  character containing the name of samples you want to modify.

- new_names:

  character containing the new name of samples.

- old_new_names_table:

  alternatively, you can bass here a table with old names in the first
  column and new names in the second column. If you pass an argument
  here, this will be considered instead of old_names and new_names.

## Value

the GCPlist with the sample names updated in the sampleINFO, and
intensities data frames.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist01 <- GCP_ChangeSampleNames(GCPlist = GCPlist00,
                                   old_names = GCPlist00$sampleINFO$Sample,
                                   new_names = c("S_1", "S_2", "S_3", "S_4", "S_5", "V_1", "V_2", "V_3", "V_4", "V_5"))

} # }

```
