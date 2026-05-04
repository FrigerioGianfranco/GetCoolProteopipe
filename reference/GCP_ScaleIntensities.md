# Scale intensities

In the intensities data table, it scale all the intensities as
specified.

## Usage

``` r
GCP_ScaleIntensities(
  GCPlist,
  subtract = c("shift_median", "mean", "median"),
  divide = c("sqrt_sd", "sd", "maxmin"),
  by_sample = TRUE,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups")
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- subtract:

  NULL, or "shift_median", or "mean", or "median". Specify it to
  subtract the mean or median value to all values. If "shift_median", an
  offset will be subtracted in order to center the medians to the global
  median.

- divide:

  NULL, or "sqrt_sd", or "sd", or "maxmin". Specify it to divide all
  values by the standard deviation, by the square root of the standard
  deviation, or by the range (max-min).

- by_sample:

  logical. If TRUE, the centering and scaling will be column-wise (per
  sample); if FALSE; row-wise (per protein).

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups. It needs to be passed
  only if subtract is "shift_median" and you if want to center each
  sample in a group to the global median of that group, instead of the
  overall global median.

## Value

a GCPlist list with the scaled data intensity tables.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist10 <- GCP_ScaleIntensities(GCPlist = GCPlist09,
                                  subtract = "shift_median",
                                  divide = "sqrt_sd",
                                  by_sample = TRUE)


} # }

```
