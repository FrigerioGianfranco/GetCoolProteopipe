# Get Principal Component analysis.

It performs a principal component analysis.

## Usage

``` r
GCP_PCA(GCPlist, center = TRUE, scale. = FALSE)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- center:

  logical. Whether the variables should be shifted to be zero centered
  (as in the prcomp function).

- scale.:

  logical. whether the variables should be scaled to have unit variance
  before the analysis takes place (as in prcomp function).

## Value

The GCPlist will be returned with the scores in the sampleINFO and the
loadings in the proteinINFO.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist12 <- GCP_PCA(GCPlist = GCPlist11)

} # }


```
