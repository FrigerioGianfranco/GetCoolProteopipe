# Remove samples from the dataset

Remove samples from both the intensities and SampleINFO data frames.
Only one of the arguments among remove_samples, keep_samples,
remove_groups, or keep_groups has to be specified.

## Usage

``` r
GCP_RemoveSamples(
  GCPlist,
  remove_samples = NULL,
  keep_samples = NULL,
  remove_groups = NULL,
  keep_groups = NULL,
  name_column_groups = NULL
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- remove_samples:

  NULL or a character containing the name of samples to remove.

- keep_samples:

  NULL or a character containing the name of samples to keep.

- remove_groups:

  NULL or a character containing the name of groups to remove.

- keep_groups:

  NULL or a character containing the name of groups to keep.

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups. It must be specified if
  you used remove_groups or keep_groups.

## Value

the GCPlist with samples removed from the sampleINFO, and intensities
data frames.

## Examples

``` r
if (FALSE) { # \dontrun{

# specify what to remove:

GCPlist00rm1 <- GCP_RemoveSamples(GCPlist = GCPlist00,
                                  remove_samples = c("S2", "S5"),
                                  keep_samples = NULL)

# or specify what to keep (same output as above):

GCPlist00rm2 <- GCP_RemoveSamples(GCPlist = GCPlist00,
                                  remove_samples = NULL,
                                  keep_samples = c("S1", "S3", "S4", "V1", "V2", "V3", "V4", "V5"))

} # }

```
