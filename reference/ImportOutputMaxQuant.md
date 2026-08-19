# Importing the output of MaxQuant

Starting from the output of the MaxQuant, it performs some cleaning such
as filtering out proteins that are only identified by site, the reverse
and the potential contaminants, then it adds protein names from a fasta
database; finally it creates the GCP list that will be used throughout
this pipeline.

## Usage

``` r
ImportOutputMaxQuant(
  MaxQuant_table_name,
  samples_info = NULL,
  fasta_database = NA,
  raw_or_LFQ = c("LFQ", "raw"),
  prioritize_MaxQuant_names = TRUE,
  remove_identified_by_site = TRUE,
  remove_reverse = TRUE,
  remove_potential_contaminant = TRUE
)
```

## Arguments

- MaxQuant_table_name:

  a character vector of length 1 with the name of the MaxQuant table
  file in the current working directory, which must be in the .txt
  format.

- samples_info:

  NULL or NA or a character vector of length 1 with the name of the
  table in the current working directory, containing information for
  each sample. The table must be in txt, csv, or xslsx format. In
  particular, the first column of the table must contain the names of
  the samples exactly as they are in the MaxQuant table.

- fasta_database:

  NULL or NA or either "human" or "mouse", or a name of a table in the
  current working directory (in .txt or .csv format). You can specify
  the fasta database to use to fill the missing Protein names and Gene
  names from the MaxQuant table. If you choose "mouse" or "human", the
  fasta table implemented were downloaded and reprocessed from Mascot on
  5 May 2024.

- raw_or_LFQ:

  one of the following: "raw", "LFQ". Only such data intensities will be
  imported.

- prioritize_MaxQuant_names:

  logical. If TRUE and if a fasta_database is provided, the final
  "Protein names" and "Gene names" will be primarily taken from the
  MaxQuant table (they will be taken from the fasta database only if
  missing). If FALSE, the opposite will happen.

- remove_identified_by_site:

  logical. Do you want to remove rows that contains "+" in the column
  "Only identified by site"?

- remove_reverse:

  logical. Do you want to remove rows that contains "+" in the column
  "Reverse"?

- remove_potential_contaminant:

  logical. Do you want to remove rows that contains "+" in the column
  "Potential contaminant"?

## Value

a GCPlist, i.e.: a list with 3 dataframes (tibbles):

- `intensities`: the protein intensities.

- `proteinINFO`: all the information for each rows.

- `sampleINFO`: the information about the samples.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist00 <- ImportOutputMaxQuant(MaxQuant_table_name = "MaxQuant OUTPUT FILE NAME.txt",     ## you can just put the name of your file in your current working directory
                                  samples_info = NULL,                                       ## also here potentially
                                  fasta_database = "mouse",
                                  raw_or_LFQ = "raw",
                                  prioritize_MaxQuant_names = TRUE,
                                  remove_identified_by_site = TRUE,
                                  remove_reverse = TRUE,
                                  remove_potential_contaminant = TRUE)


} # }

```
