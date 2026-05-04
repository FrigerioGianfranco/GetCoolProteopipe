# Log-transform intensities

It log-transforms all the intensities of the intensities table.

## Usage

``` r
GCP_LogTransformIntensities(
  GCPlist,
  base = getOption("GetCoolProteopipe.log_base")
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- base:

  numeric. The base of the logarithm to use.

## Value

a GCPlist list with log-transformed data intensity table.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist06 <- GCP_LogTransformIntensities(GCPlist = GCPlist05,
                                         base = 2)

} # }

```
