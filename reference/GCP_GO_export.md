# Export Gene Ontology (GO) results to Excel.

It exports the results obtained with the function GCP_GO() in one or two
.xlsx files in the current working directory.

## Usage

``` r
GCP_GO_export(
  GO_list,
  file_name_GO = "GO_results",
  export_considered_gene_names = FALSE,
  file_name_considered_gene_names = "GO_considered_gene_names"
)
```

## Arguments

- GO_list:

  A list obtained with the function GCP_GO().

- file_name_GO:

  The name you want to assign to the generated file, or the full file
  path.

- export_considered_gene_names:

  logical. If TRUE, an additional file will be exported, with the
  information related to genes used in the computed GOs.

- file_name_considered_gene_names:

  The name you want to assign to the file related to the considered gene
  names, or the full file path, if export_considered_gene_names is TRUE.

## Value

It exports one or two file in the current working directory.

## Examples

``` r

if (FALSE) { # \dontrun{

this_GO_list <- GCP_GO(GCPlist = GCPlist12a2,
                       name_column_groups = "Group_multi",
                       OrgDb = "mouse")

GCP_GO_export(GO_list = this_GO_list,
              file_name_GO = "GO_results",
              export_considered_gene_names = TRUE,
              file_name_considered_gene_names = "GO_considered_gene_names")

} # }
```
