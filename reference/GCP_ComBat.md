# Adjust for batch effects

It performs adjustment for batch effects using an empirical Bayes
framework applying the function ComBat from the sva package.

## Usage

``` r
GCP_ComBat(GCPlist, batch = "Batch", ...)
```

## Arguments

- GCPlist:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- batch:

  character of length 1 OR a numeric/factor vector. The name of the
  column of the sampleINFO table containing the batch indications OR a
  numeric or factor vector containing the batch indications, that will
  also be added as column in the sampleINFO data frame under a new
  'batch' column.

- ...:

  Additional arguments passed to sva::ComBat.

## Value

The GCPlist with the desired intensity values adjusted for batch
effects.

## Examples

``` r
if (FALSE) { # \dontrun{


# passing a numeric/factor vector directly to the 'batch' argument:

GCPlist14b1 <- GCP_ComBat(GCPlist = GCPlist14b,
                          batch = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1))

# passing a column name of sampleINFO, containing the numeric/factor vector, to the 'batch' argument:

GCPlist14b$sampleINFO[, "batches"] <- factor(c("batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1"), levels = c("batch1", "batch2", "batch3"))

GCPlist14b2 <- GCP_ComBat(GCPlist = GCPlist14b,
                          batch = "batches")

} # }

```
