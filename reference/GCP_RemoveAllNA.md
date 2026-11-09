# Remove samples with all missing values intensities

Starting from a GCPlist, it removes proteins which intensity is missing
in all the samples.

## Usage

``` r
GCP_RemoveAllNA(GCPlist)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

## Value

the GCPlist in which each table has potentially a reduced number of
rows.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist03 <- GCP_RemoveAllNA(GCPlist02)

} # }

```
