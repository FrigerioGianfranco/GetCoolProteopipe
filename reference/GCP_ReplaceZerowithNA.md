# Replacing zero with NA

In the intensities table, it replaces all zeros with missing values.

## Usage

``` r
GCP_ReplaceZerowithNA(GCPlist)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

## Value

a GCPlist list with zero replaced with NAs in the data intensity table.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist04 <- GCP_ReplaceZerowithNA(GCPlist03)

} # }

```
