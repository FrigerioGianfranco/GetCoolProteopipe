# ratio PTM on protein intensities

It takes each valid feature intensity from PTM and it divides it by the
relative intensity of the protein; the match is based on the Accession.

## Usage

``` r
GCP_ratioPTMprotein(GCPlistPTM, GCPlistProteins, are_log_transf = TRUE)
```

## Arguments

- GCPlistPTM:

  a list initially created with ImportOutputMaxQuan(), ImportPTMs(), or
  ImportOutputProtDiscov().

- GCPlistProteins:

  a list created with the ImportOutputMaxQuant function.

- are_log_transf:

  logical. If the intensities are log-transformed, specify here as TRUE,
  so the subtraction will be performed instead of the ratio.

## Value

a GCPlist list with the calculated ratios of intensities in the
intensities table.

## Examples

``` r
if (FALSE) { # \dontrun{


GCPlistPTM09 <- GCP_ratioPTMprotein(GCPlistPTM = GCPlistPTM08,
                                    GCPlistProteins = GCPlist14,
                                    are_log_transf = TRUE)


} # }

```
