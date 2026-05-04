# Exporting the GCPlist

It exports a GPClist in a single table file in the current working
directory.

## Usage

``` r
ExportGCPlist(
  GCPlist,
  filename = "GCP.xlsx",
  exportype = c("intensities", "proteinINFO", "sampleINFO", "all", "sheets"),
  intensity_indication = TRUE,
  protgenenames = TRUE,
  protIDclean = TRUE,
  specific_columns = "all_stat"
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function or a list of
  such lists. If the second is passed, the following filename argument
  must end with ".xlsx" as a sheet for each of those list will be
  created.

- filename:

  character. The name of the file to create. It must end with ".txt",
  ".csv", or ".xlsx"; and the file will be accordingly created of that
  format.

- exportype:

  exportype one of the following: "intensities", "proteinINFO",
  "sampleINFO", "all", "sheets". If "intensities", a table with the
  intensities will be exported; if "proteinINFO", a table with all the
  proteinINFO; if "sampleINFO", a table with all the sampleINFO; if
  "all", a table with combined intensities and proteinINFO will be
  exported; if "sheets", an Excel table with 3 sheets will be exported.
  Please note that "all" might be too big to be suitably exported as a
  '.xlsx' file: if so, export it as '.txt' or '.csv'.

- intensity_indication:

  logical. if TRUE and if exportype is "intensities", the column names
  of samples will start with "Intensity ". Please note that this will
  happen anyway if exportype is "all".

- protgenenames:

  logical. If TRUE and if exportype is "intensities" or "proteinINFO",
  it will export the columns 'Protein IDs' (or 'Accession', see below),
  'Protein names', and 'Gene names'.

- protIDclean:

  logical. If TRUE (and if protgenenames is TRUE, and if exportype is
  "intensities" or "proteinINFO"), it exports the column 'Accession',
  which contains only the first code of each protein of the Protein IDs.
  If FALSE it exports the complete 'Protein IDs' column.

- specific_columns:

  character. If exportype is "intensities" or "proteinINFO", you can
  specify here some additional columns to export from the proteinINFO
  table. Moreover, you can simply pass here a single element to export
  all the related result columns, in particular: "shap_test", "PC",
  "ttest", "FC", "ANOVA", "presence_group"; and also "all_stat" for all
  of those.

## Value

Export the file in the current working directory.

## Examples

``` r
if (FALSE) { # \dontrun{

# exporting a single GCPlist:

ExportGCPlist(GCPlist14, "GCPlist14_single_export.txt")


# exporting multiple GCPlists:

GCPlist_of_lists <- list(asimported = GCPlist00,
                         `before processing` = GCPlist05,
                         `post processing` = GCPlist11,
                         with_statistics = GCPlist14,
                         with_mock_ANOVA = GCPlist14a2,
                         with_mock_batchcorr = GCPlist14b2)

ExportGCPlist(GCPlist_of_lists, "GCPlists_multiple_export.xlsx")


} # }

```
