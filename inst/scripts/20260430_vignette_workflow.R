
## Installation (only for the first time)

if (!require("devtools", quietly = TRUE)) {install.packages("devtools")}

devtools::install_github("FrigerioGianfranco/GetCoolProteopipe", dependencies = TRUE)


## Loading the package

library(GetCoolProteopipe)


# Untargeted proteomics

## Import of MaxQuant output

proteinGroups_path <- system.file("extdata", "proteinGroups_vignettes.txt.gz", package = "GetCoolProteopipe")

GCPlist00 <- ImportOutputMaxQuant(MaxQuant_table_name = proteinGroups_path,   ## you can just put the name of your file in your current working directory
                                  samples_info = NULL,                       ## also here potentially
                                  fasta_database = "mouse",
                                  raw_or_LFQ = "raw",
                                  prioritize_MaxQuant_names = TRUE,
                                  remove_identified_by_site = TRUE,
                                  remove_reverse = TRUE,
                                  remove_potential_contaminant = TRUE)

## Cleaning of imported data

GCPlist01 <- GCP_RemoveAllZero(GCPlist00)


GCPlist02 <- GCP_ReplaceZerowithNA(GCPlist01)


## Managing samples


### GCP_RemoveSamples, specify what to remove:

GCPlist02rm1 <- GCP_RemoveSamples(GCPlist = GCPlist02,
                                  remove_samples = c("S2", "S5"),
                                  keep_samples = NULL)

### GCP_RemoveSamples, or specify what to keep (same output as above):

GCPlist02rm2 <- GCP_RemoveSamples(GCPlist = GCPlist02,
                                  remove_samples = NULL,
                                  keep_samples = c("S1", "S3", "S4", "V1", "V2", "V3", "V4", "V5"))



GCPlist03 <- GCP_ChangeSampleNames(GCPlist = GCPlist02,
                                   old_names = GCPlist02$sampleINFO$Sample,
                                   new_names = c("S_1", "S_2", "S_3", "S_4", "S_5", "V_1", "V_2", "V_3", "V_4", "V_5"))

GCPlist04 <- GCP_AssignGroups(GCPlist = GCPlist03,
                              automatic_assignment = "groupfirst",
                              separator_automatic_assignment = "_",
                              name_column_groups = "Condition",
                              controlgroup = "S")

set_name_column_groups("Condition")

set_col_pal(c(S = "green", V = "orange"))


### GCP_ReorderSamples, based on the groups (by default based on the set name_column_groups):

GCPlist04rd1 <- GCP_ReorderSamples(GCPlist04)

### GCP_ReorderSamples, specifying the samples:

GCPlist04rd2 <- GCP_ReorderSamples(GCPlist04, sample_names_ordered = c("S_1", "V_1", "S_2", "V_2", "S_3", "V_3", "S_4", "V_4", "S_5", "V_5"))




## Processing

GCPlist05 <- GCP_LogTransformIntensities(GCPlist = GCPlist04)


GCPlist06 <- GCP_FilterNAperCondition(GCPlist = GCPlist05,
                                      ratio = 0.5)



GCPlist07 <- GCP_ProteinsGrouped(GCPlist06)


GCPlist08 <- GCP_NAimputation(GCPlist = GCPlist07,
                              quant_rate = 0.5)



GCPlist09 <- GCP_ScaleIntensities(GCPlist = GCPlist08,
                                  subtract = "shift_median",
                                  divide = "sqrt_sd",
                                  by_sample = TRUE)


### Batch correction:

GCPlist09b <- GCPlist09

#### passing a numeric/factor vector directly to the 'batch' argument:

GCPlist09b1 <- GCP_ComBat(GCPlist = GCPlist09b,
                          batch = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1))

#### passing a column name of sampleINFO, containing the numeric/factor vector, to the 'batch' argument:

GCPlist09b$sampleINFO[, "batches"] <- factor(c("batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1"), levels = c("batch1", "batch2", "batch3"))

GCPlist09b2 <- GCP_ComBat(GCPlist = GCPlist09b,
                          batch = "batches")



## Visualisation of data before and after processing



## before processing:
GCPlist04tn <- GCP_TestNormality(GCPlist04)

## after processing:
GCPlist09tn <- GCP_TestNormality(GCPlist09)



## before processing:
Fig01_IntensityDistribution_before_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist04,
                                                                            Title = "Distribution of intensities, before processing")
export_figures(Fig01_IntensityDistribution_before_processing)    ## this would create a .png file with that name on your current working directory

print(Fig01_IntensityDistribution_before_processing)

## after processing:
Fig02_IntensityDistribution_post_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist09,
                                                                          Title = "Distribution of intensities, post processing")
export_figures(Fig02_IntensityDistribution_post_processing)

print(Fig02_IntensityDistribution_post_processing)



## before processing:
Fig03_QQ_plot_before_processing <- GCP_QQplotIntensities(GCPlist = GCPlist04,
                                                         Title = "QQ plot - before processing")
export_figures(Fig03_QQ_plot_before_processing)

print(Fig03_QQ_plot_before_processing)

## after processing:
Fig04_QQ_plot_post_processing <- GCP_QQplotIntensities(GCPlist = GCPlist09,
                                                       Title = "QQ plot - post processing")
export_figures(Fig04_QQ_plot_post_processing)

print(Fig04_QQ_plot_post_processing)




Fig05_BarPlot_post_filteringNA <- GCP_BarPlot(GCPlist = GCPlist07,
                                              label_numbers = TRUE,
                                              showCV = TRUE,
                                              rotate_sample_names = TRUE)
export_figures(Fig05_BarPlot_post_filteringNA)

print(Fig05_BarPlot_post_filteringNA)

Fig06_Venn_post_filteringNA <- GCP_Venn(GCPlist07)
export_figures(Fig06_Venn_post_filteringNA)

print(Fig06_Venn_post_filteringNA)



SummaryTable_ProteinsGrouped_filteringNA <- GCP_ProteinsGroupedSummary(GCPlist07)
export_the_table(SummaryTable_ProteinsGrouped_filteringNA)   ## this would quickly export that table in your working directory!

print(SummaryTable_ProteinsGrouped_filteringNA)






## before processing:
Fig07_BoxPlot_before_processing <- GCP_BoxPlots(GCPlist = GCPlist04,
                                                by_samples = TRUE,
                                                Title = "Distribution of intensities, before processing")
export_figures(Fig07_BoxPlot_before_processing)

print(Fig07_BoxPlot_before_processing)

## after processing:
Fig08_BoxPlot_post_processing <- GCP_BoxPlots(GCPlist = GCPlist09,
                                              by_samples = TRUE,
                                              Title = "Distribution of intensities, post processing")
export_figures(Fig08_BoxPlot_post_processing)

print(Fig08_BoxPlot_post_processing)



## before processing:
Fig09_Dendogram_before_processing <- GCP_Dendogram(GCPlist04)
export_figures(Fig09_Dendogram_before_processing)

print(Fig09_Dendogram_before_processing)


## after processing:
Fig10_Dendogram_post_processing <- GCP_Dendogram(GCPlist09)
export_figures(Fig10_Dendogram_post_processing)

print(Fig10_Dendogram_post_processing)





Fig11_PCA_score_plot <- GCP_plotPCA(GCPlist = GCPlist09,
                                    scores_or_loadings = "scores",
                                    PC_to_plot = c("PC1", "PC2"),
                                    name_column_labels = "Sample",
                                    ellipses_on_score = TRUE)
export_figures(Fig11_PCA_score_plot)

print(Fig11_PCA_score_plot)

Fig12_PCA_loading_plot <- GCP_plotPCA(GCPlist = GCPlist09,
                                      scores_or_loadings = "loadings",
                                      PC_to_plot = c("PC1", "PC2"),
                                      name_column_groups_loading = NULL,
                                      name_column_labels_loading = "Protein names",
                                      col_pal_loading = NULL,
                                      ellipses_on_loading = FALSE)
export_figures(Fig12_PCA_loading_plot)

print(Fig12_PCA_loading_plot)


GCPlist10 <- GCP_PCA(GCPlist = GCPlist09)





## Statistics

GCPlist11 <- GCP_ttest(GCPlist = GCPlist10,
                       paired = TRUE,
                       FDR = TRUE,
                       pcutoff = 0.05)


GCPlist12 <- GCP_FoldChange(GCPlist = GCPlist11,
                            paired = TRUE,
                            are_log_transf = TRUE)




GCPlist12a <- GCPlist12
GCPlist12a$sampleINFO[, "Group_multi"] <- factor(c("gr1", "gr1", "gr2", "gr2", "gr3", "gr3", "gr3", "gr4", "gr4", "gr4"), levels = c("gr1", "gr2", "gr3", "gr4"))


GCPlist12a1 <- GCP_ANOVA(GCPlist = GCPlist12a,
                         name_column_groups = "Group_multi")


GCPlist12a2 <- GCP_FoldChange_Multi(GCPlist = GCPlist12a1,
                                    name_column_groups = "Group_multi",
                                    paired = FALSE,
                                    are_log_transf = TRUE)




## Visualisation of singificant proteins


Fig13_Volcano_plot <- GCP_Volcano(GCPlist = GCPlist12,
                                  x_val = "logFC",
                                  y_val = "ttest_Pvalues",
                                  pcutoff_colored = 0.05,
                                  pcutoff_line = 0.05,
                                  pcutoff_prot_label = 0.005,
                                  name_column_proteinlabels = "Protein names",
                                  name_column_proteingroups = NULL,
                                  col_pal_difference = c("grey", "blue", "red"),
                                  col_pal_groups = NULL)
export_figures(Fig13_Volcano_plot)

print(Fig13_Volcano_plot)







GCPlist12f <- GCP_FilterProteins(GCPlist = GCPlist12,
                                 operation = "ttest_Pvalues < 0.05")


Fig14_The_Heat_Map_sign <- GCP_HeatMap(GCPlist = GCPlist12f,
                                       name_column_groups = c("Condition", "Replicate"),
                                       name_column_labels = "Sample",
                                       name_column_groups_protein = NULL,
                                       name_column_labels_protein = "Protein names",
                                       name_rows = FALSE,
                                       name_columns = TRUE,
                                       rotate_name_columns = TRUE,
                                       col_pal_list = list(Condition = c(S = "green", V = "blue"),
                                                           Replicate = c("lightskyblue", "lightskyblue1", "lightskyblue2", "lightskyblue3", "lightskyblue4")))
export_figures(Fig14_The_Heat_Map_sign)

print(Fig14_The_Heat_Map_sign)



Fig15_BoxPlots_sign <- GCP_BoxPlots(GCPlist = GCPlist12f,
                                    by_samples = FALSE)
export_figures(Fig15_BoxPlots_sign)                            ## exporting each box-plot as a single png file
export_figures(Fig15_BoxPlots_sign, exprt_fig_type = "pdf")    ## exporting all box-plots in a pdf file with many pages

print(Fig15_BoxPlots_sign[1:4])



## export

## exporting a single GCPlist:

ExportGCPlist(GCPlist12, "GCPlist12_single_export.txt")


## exporting multiple GCPlists:

GCPlist_of_lists <- list(asimported = GCPlist00,
                         `before processing` = GCPlist04,
                         `post processing` = GCPlist09,
                         with_statistics = GCPlist12,
                         with_mock_ANOVA = GCPlist12a2,
                         with_mock_batchcorr = GCPlist09b2)

ExportGCPlist(GCPlist_of_lists, "GCPlists_multiple_export.xlsx")








# Post-translational modifications (PTMs)


proteinGroupsPTM_path <- system.file("extdata", "AcetylKSites_vignettes.txt.gz", package = "GetCoolProteopipe")

GCPlistPTM00 <- ImportPTMs(MaxQuant_table_name = proteinGroupsPTM_path,     ## you can just put the name of your file in your current working directory
                           samples_info = NULL,                             ## also here potentially
                           fasta_database = "mouse",
                           prioritize_MaxQuant_names = TRUE,
                           pattern_intensity = "Intensity ",
                           pattern_isoforms = "___")

GCPlistPTM01 <- GCP_ChangeSampleNames(GCPlist = GCPlistPTM00,
                                      old_names = GCPlistPTM00$sampleINFO$Sample,
                                      new_names = c("S_1", "S_2", "S_3", "S_4", "S_5", "V_1", "V_2", "V_3", "V_4", "V_5"))


GCPlistPTM02 <- GCP_AssignGroups(GCPlist = GCPlistPTM01,
                                 automatic_assignment = "groupfirst",
                                 separator_automatic_assignment = "_",
                                 name_column_groups = "Condition",
                                 controlgroup = "S")

set_name_column_groups("Condition")


GCPlistPTM03 <- GCP_RemoveAllZero(GCPlistPTM02)

GCPlistPTM04 <- GCP_ReplaceZerowithNA(GCPlistPTM03)

GCPlistPTM05 <- GCP_LogTransformIntensities(GCPlist = GCPlistPTM04)

GCPlistPTM06 <- GCP_FilterNAperCondition(GCPlist = GCPlistPTM05,
                                         ratio = 0.5)

GCPlistPTM07 <- GCP_NAimputation(GCPlist = GCPlistPTM06,
                                 quant_rate = 0.5)


GCPlistPTM08 <- GCP_ScaleIntensities(GCPlist = GCPlistPTM07,
                                     subtract = "shift_median",
                                     divide = "sqrt_sd",
                                     by_sample = TRUE)

GCPlistPTM09 <- GCP_ratioPTMprotein(GCPlistPTM = GCPlistPTM08,
                                    GCPlistProteins = GCPlist12,
                                    are_log_transf = TRUE)

GCPlist_of_PTMlists <- list(asimported = GCPlistPTM00,
                            `before processing` = GCPlistPTM04,
                            `post processing` = GCPlistPTM08,
                            with_ratio_protein = GCPlistPTM09)
ExportGCPlist(GCPlist_of_PTMlists, "GCPlistsPTM_multiple_export.xlsx")

