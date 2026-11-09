# Performing Gene Ontology (GO) Enrichment Analyses

It performs multiple Gene Ontology (GO) Enrichment Analyses for each
pairwise comparison using the function enrichGO() from the
clusterProfiler package. The genes used are the ones retrieved from the
'Gene names' column of proteinINFO.

## Usage

``` r
GCP_GO(
  GCPlist,
  name_column_groups = getOption("GetCoolProteopipe.name_column_groups"),
  use_FDRPvalue = TRUE,
  Pvalue_column = NULL,
  Pvalue_cutoff = 0.05,
  differential_analysis = TRUE,
  use_logFC = TRUE,
  FC_column = NULL,
  FC_cutoff = 0,
  background_from_data = TRUE,
  OrgDb = "human",
  keyType = "SYMBOL",
  ont = c("ALL", "BP", "MF", "CC"),
  pvalueCutoff = 0.05,
  pAdjustMethod = "none",
  qvalueCutoff = 0.2,
  ...
)
```

## Arguments

- GCPlist:

  a list created with the ImportOutputMaxQuant function.

- name_column_groups:

  NULL or character of length 1. The name of the column of the
  sampleINFO table containing the sample groups. If NULL, you will need
  to specify Pvalue_column and, possibly, FC_column, see below.

- use_FDRPvalue:

  logical. If name_column_groups is specified and if this is TRUE, this
  function will consider columns with FDR p-values in proteinINFO;
  otherwise, it will consider columns with the p-values in proteinINFO.

- Pvalue_column:

  NULL or character. If name_column_groups was not specified, you can
  directly pass here the column names of proteinINFO you want to
  consider as p-values.

- Pvalue_cutoff:

  numeric. P-values below this numner will be consider as statistically
  significant and relative gene names will be considered.

- differential_analysis:

  logical. If TRUE, it will perform two separate analyses for each
  pairwise condition considering the Fold Changes, see below;if FALSE;
  it will only perform an analyses for each pairwise condition.

- use_logFC:

  logical. If name_column_groups is specified and if this is TRUE, this
  function will consider columns with log FC in proteinINFO; otherwise,
  it will consider columns with FC in proteinINFO.

- FC_cutoff:

  numeric. significant FC higher than this value will be considered for
  the differential 'up' analysis, and significant FC lower than this
  value will be considered for the differential 'down' analysis.

- background_from_data:

  logical. If TRUE; the "universe" passed to enrichGO will be valid
  unique gene names from all the not-NA protein groups; if FALSE, all
  genes listed in the database will be used as background.

- OrgDb:

  One of "human", "mouse", or "rat"; or an OrgDb object to pass to
  clusterProfiler::enrichGO. If "human", the OrgDb "org.Hs.eg.db" will
  be used; if "mouse", "org.Mm.eg.db"; if "rat", "org.Rn.eg.db". In any
  case, you need to have that package installed; if you don't, please
  follow the instruction explained in the error message to do so.

- keyType:

  character of length 1. The keyType to pass to
  clusterProfiler::enrichGO. Considering that the values are taken from
  the 'Gene names' of proteinINFO, this should be "SYMBOL" as the
  defalut.

- ont:

  One of "ALL", "BP", "MF", or "CC". ont to pass to
  clusterProfiler::enrichGO.

- pvalueCutoff:

  numeric. pvalueCutoff to pass to clusterProfiler::enrichGO.

- pAdjustMethod:

  one of "none", "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  "fdr". pAdjustMethod to pass to clusterProfiler::enrichGO.

- qvalueCutoff:

  numeric. qvalueCutoff to pass to clusterProfiler::enrichGO.

- ...:

  Further argument to pass to clusterProfiler::enrichGO.

## Value

a list with three elements:

- `considered_gene_names`: a list with all the gene names considered.

- `SUMMARY_gene_names`: a tibble which is a summary of the numbers of
  gene names considered.

- `GO_results`: A list with dataframes, each containing the results of
  the GO.

## Examples

``` r
if (FALSE) { # \dontrun{



# specifying a name_column_groups:

the_list_of_GO_results <- GCP_GO(GCPlist = GCPlist12,
                                 name_column_groups = "Condition",
                                 use_FDRPvalue = TRUE,
                                 Pvalue_cutoff = 0.05,
                                 differential_analysis = TRUE,
                                 use_logFC = TRUE,
                                 FC_cutoff = 0,
                                 background_from_data = TRUE,
                                 OrgDb = "mouse")

# not specifying a name_column_groups:

the_list_of_GO_results <- GCP_GO(GCPlist = GCPlist12,
                                 name_column_groups = NULL,
                                 Pvalue_column = c("ttest_V_vs_S_PvaluesFDR"),
                                 Pvalue_cutoff = 0.05,
                                 differential_analysis = TRUE,
                                 FC_column = c("V_vs_S_logFC"),
                                 FC_cutoff = 0,
                                 background_from_data = TRUE,
                                 OrgDb = "mouse")


} # }
```
