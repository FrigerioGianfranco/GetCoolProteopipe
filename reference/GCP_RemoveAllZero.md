# Remove samples with all zero intensities

Starting from a GCPlist, it removes proteins which intensity is equal to
zero in all the samples.

## Usage

``` r
GCP_RemoveAllZero(GCPlist)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

## Value

the GCPlist in which each table has potentially a reduced number of
rows.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist03 <- GCP_RemoveAllZero(GCPlist02)

} # }

```
