# Importing the output of Proteome Discoverer

Starting from the output of the Proteome Discoverer, it performs some
cleaning such as removing samples with all missing values and/or
filtering out proteins that are contaminants, then it can also add
protein names from a fasta database; finally it creates the GCP list
that will be used throughout this pipeline.

## Usage

``` r
ImportOutputProtDiscov(
  ProtDiscov_table_name,
  ProtDiscov_InputFiles = NULL,
  samples_info = NULL,
  raw_or_norm = c("raw", "norm"),
  restore_sample_names = TRUE,
  remove_empty_columns = TRUE,
  fasta_database = NA,
  prioritize_ProtDiscov_names = TRUE,
  remove_contaminants = TRUE
)
```

## Arguments

- ProtDiscov_table_name:

  a character vector of length 1 with the name of the Proteome
  Discoverer table file, as exported, in the current working directory,
  which must be in the .txt format.

- ProtDiscov_InputFiles:

  NULL or a character vector of length 1 with the name of the Proteome
  Discoverer input file table, as exported, in the current working
  directory, which must be in the .txt format.

- samples_info:

  NULL or NA or a character vector of length 1 with the name of the
  table in the current working directory, containing information for
  each sample. The table must be in txt, csv, or xslsx format. In
  particular, the first column of the table must contain the names of
  the samples exactly the same considered.

- raw_or_norm:

  either "raw" or "norm". Choose whether you wish to import the raw
  abundances or the normalised ones.

- restore_sample_names:

  logical. If TRUE, you must provide the ProtDiscov_InputFiles and the
  original sample names will be considered. If FALSE, the Proteome
  Discoverer sample names (like F1, F2, ...) will be considered.

- remove_empty_columns:

  logical. If TRUE, samples that have full missing values abundances
  will be removed.

- fasta_database:

  NULL or NA or either "human" or "mouse", or a name of a table in the
  current working directory (in .txt or .csv format). You can specify
  the fasta database to use to fill the missing Protein names and Gene
  names from the Proteome Discoverer table. If you choose "mouse" or
  "human", the fasta table implemented were downloaded and reprocessed
  from Mascot on 5 May 2024.

- prioritize_ProtDiscov_names:

  logical. If TRUE and if a fasta_database is provided, the final
  "Protein names" and "Gene names" will be primarily taken from the
  Proteome Discoverer table (they will be taken from the fasta database
  only if missing). If FALSE, the opposite will happen.

- remove_contaminants:

  logical. Do you want to remove rows that are TRUE in the column
  "Contaminant"?

## Value

a GCPlist, i.e.: a list with 3 dataframes (tibbles):

- `intensities`: the protein intensities.

- `proteinINFO`: all the information for each rows.

- `sampleINFO`: the information about the samples.

## Examples

``` r
if (FALSE) { # \dontrun{

GCPlist00P <- ImportOutputProtDiscov(ProtDiscov_table_name = "Proteome Discoverer OUTPUT FILE NAME.txt",     ## put here the name of the Proteome Discoverer output table present in your current working directory
                                     ProtDiscov_InputFiles = "Proteome Discoverer INPUT FILE NAME.txt",      ## you could also put the name the input file if you want to restore the original sample names
                                     samples_info = NULL,                                                 ## you could put here the file name of a table with further information about your samples
                                     raw_or_norm = "raw",
                                     restore_sample_names = TRUE,
                                     remove_empty_columns = TRUE,
                                     fasta_database = "human",
                                     prioritize_ProtDiscov_names = TRUE,
                                     remove_contaminants = TRUE)

} # }
```
