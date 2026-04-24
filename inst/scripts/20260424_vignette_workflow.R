
# package installation (only for the first time)

if (!require("devtools", quietly = TRUE)) {install.packages("devtools")}

devtools::install_github("FrigerioGianfranco/GetCoolProteopipe", dependencies = TRUE)


##
#### Vignette full workflow
##

library(GetCoolProteopipe)


proteinGroups_path <- system.file("extdata", "proteinGroups_vignettes.txt.gz", package = "GetCoolProteopipe")

GCPlist00 <- ImportOutputMaxQuant(MaxQuant_table_name = proteinGroups_path,   ## you can just put the name of your file in your current working directory
                                  samples_info = NULL,                       ## also here potentially
                                  fasta_database = "mouse",
                                  raw_or_LFQ = "raw",
                                  prioritize_MaxQuant_names = TRUE,
                                  remove_identified_by_site = TRUE,
                                  remove_reverse = TRUE,
                                  remove_potential_contaminant = TRUE)



# improving names / groups

## GCP_RemoveSamples, specify what to remove:

GCPlist00rm1 <- GCP_RemoveSamples(GCPlist = GCPlist00,
                                  remove_samples = c("S2", "S5"),
                                  keep_samples = NULL)

## GCP_RemoveSamples, or specify what to keep (same output as above):

GCPlist00rm2 <- GCP_RemoveSamples(GCPlist = GCPlist00,
                                  remove_samples = NULL,
                                  keep_samples = c("S1", "S3", "S4", "V1", "V2", "V3", "V4", "V5"))



GCPlist01 <- GCP_ChangeSampleNames(GCPlist = GCPlist00,
                                   old_names = GCPlist00$sampleINFO$Sample,
                                   new_names = c("S_1", "S_2", "S_3", "S_4", "S_5", "V_1", "V_2", "V_3", "V_4", "V_5"))

GCPlist02 <- GCP_AssignGroups(GCPlist = GCPlist01,
                              automatic_assignment = "groupfirst",
                              separator_automatic_assignment = "_",
                              name_column_groups = "Condition",
                              controlgroup = "S")

set_name_column_groups("Condition")

set_col_pal(c(S = "green", V = "orange"))


## GCP_ReorderSamples, based on the groups (by default based on the set name_column_groups):

GCPlist02rd1 <- GCP_ReorderSamples(GCPlist02)

## GCP_ReorderSamples, specifying the samples:

GCPlist02rd2 <- GCP_ReorderSamples(GCPlist02, sample_names_ordered = c("S_1", "V_1", "S_2", "V_2", "S_3", "V_3", "S_4", "V_4", "S_5", "V_5"))




# pre-processing

GCPlist03 <- GCP_RemoveAllZero(GCPlist02)


GCPlist04 <- GCP_ReplaceZerowithNA(GCPlist03)


# normality test and visualisation pre-processing:

GCPlist05 <- GCP_TestNormality(GCPlist04)



Fig01_IntensityDistribution_before_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist05,
                                                                            Title = "Distribution of intensities, before processing")
export_figures(Fig01_IntensityDistribution_before_processing)


Fig02_QQ_plot_before_processing <- GCP_QQplotIntensities(GCPlist = GCPlist05,
                                                         Title = "QQ plot - before processing")
export_figures(Fig02_QQ_plot_before_processing)




Fig03_BoxPlot_before_processing <- GCP_BoxPlots(GCPlist = GCPlist05,
                                                by_samples = TRUE,
                                                Title = "Distribution of intensities, before processing")
export_figures(Fig03_BoxPlot_before_processing)



Fig04_Dendogram_before_processing <- GCP_Dendogram(GCPlist05)
export_figures(Fig04_Dendogram_before_processing)



Fig05_BarPlot_before_processing <- GCP_BarPlot(GCPlist05)
export_figures(Fig05_BarPlot_before_processing)



Fig06_Venn_before_processing <- GCP_Venn(GCPlist05)
export_figures(Fig06_Venn_before_processing)




# processing:


GCPlist06 <- GCP_LogTransformIntensities(GCPlist = GCPlist05)

GCPlist07 <- GCP_FilterNAperCondition(GCPlist = GCPlist06,
                                      ratio = 0.5)



# bar plot and Venn after filteringNAperCondition:

Fig07_BarPlot_post_filteringNA <- GCP_BarPlot(GCPlist = GCPlist07,
                                              label_numbers = TRUE,
                                              showCV = TRUE,
                                              rotate_sample_names = TRUE)

export_figures(Fig07_BarPlot_post_filteringNA)

Fig08_Venn_post_filteringNA <- GCP_Venn(GCPlist07)
export_figures(Fig08_Venn_post_filteringNA)


## proteins in groups

GCPlist08 <- GCP_ProteinsGrouped(GCPlist07)

SummaryTable_ProteinsGrouped_filteringNA <- GCP_ProteinsGroupedSummary(GCPlist08)
export_the_table(SummaryTable_ProteinsGrouped_filteringNA)



# NA imputation and scaling

GCPlist09 <- GCP_NAimputation(GCPlist = GCPlist08,
                              quant_rate = 0.5)



GCPlist10 <- GCP_ScaleIntensities(GCPlist = GCPlist09,
                                  subtract = "shift_median",
                                  divide = "sqrt_sd",
                                  by_sample = TRUE)





# visualisation post processing:

GCPlist11 <- GCP_TestNormality(GCPlist10)

Fig09_IntensityDistribution_post_processing <- GCP_DensityplotIntensities(GCPlist = GCPlist11,
                                                                          Title = "Distribution of intensities, post processing")
export_figures(Fig09_IntensityDistribution_post_processing)


Fig10_QQ_plot_post_processing <- GCP_QQplotIntensities(GCPlist = GCPlist11,
                                                       Title = "QQ plot - post processing")
export_figures(Fig10_QQ_plot_post_processing)



Fig11_BoxPlot_post_processing <- GCP_BoxPlots(GCPlist = GCPlist11,
                                              by_samples = TRUE,
                                              Title = "Distribution of intensities, post processing")
export_figures(Fig11_BoxPlot_post_processing)


Fig12_Dendogram_post_processing <- GCP_Dendogram(GCPlist11)
export_figures(Fig12_Dendogram_post_processing)



# PCA on all proteins:

## to get a score plot:

Fig13_PCA_score_plot <- GCP_plotPCA(GCPlist = GCPlist11,
                                    scores_or_loadings = "scores",
                                    PC_to_plot = c("PC1", "PC2"),
                                    name_column_labels = "Sample",
                                    ellipses_on_score = TRUE)
export_figures(Fig13_PCA_score_plot)


## to get a loading plot:

Fig14_PCA_loading_plot <- GCP_plotPCA(GCPlist = GCPlist11,
                                      scores_or_loadings = "loadings",
                                      PC_to_plot = c("PC1", "PC2"),
                                      name_column_groups_loading = NULL,
                                      name_column_labels_loading = "Protein names",
                                      col_pal_loading = NULL,
                                      ellipses_on_loading = FALSE)
export_figures(Fig14_PCA_loading_plot)


## to get scores and loadings results:

GCPlist12 <- GCP_PCA(GCPlist = GCPlist11)



# t test and FC:

GCPlist13 <- GCP_ttest(GCPlist = GCPlist12,
                       paired = TRUE,
                       FDR = TRUE,
                       pcutoff = 0.05)


GCPlist14 <- GCP_FoldChange(GCPlist = GCPlist13,
                            paired = TRUE,
                            are_log_transf = TRUE)


Fig15_Volcano_plot <- GCP_Volcano(GCPlist = GCPlist14,
                                  x_val = "logFC",
                                  y_val = "ttest_PvaluesFDR",
                                  pcutoff_colored = 0.05,
                                  pcutoff_line = 0.05,
                                  pcutoff_prot_label = 0.005,
                                  name_column_proteinlabels = "Protein names",
                                  name_column_proteingroups = NULL,
                                  col_pal_difference = c("grey", "blue", "red"),
                                  col_pal_groups = NULL)
export_figures(Fig15_Volcano_plot)



# selecting only significant proteins to do the Heat Map and PCA only on significant:

GCPlist14f <- GCP_FilterProteins(GCPlist = GCPlist14,
                                 operation = "ttest_Pvalues < 0.05")



## HeatMap only on significant proteins:

Fig16_The_Heat_Map_sign <- GCP_HeatMap(GCPlist = GCPlist14f,
                                       name_column_groups = c("Condition", "Replicate"),
                                       name_column_labels = "Sample",
                                       name_column_groups_protein = NULL,
                                       name_column_labels_protein = "Protein names",
                                       name_rows = FALSE,
                                       name_columns = TRUE,
                                       rotate_name_columns = TRUE,
                                       col_pal_list = list(Condition = c(S = "green", V = "orange"),
                                                           Replicate = c("lightskyblue", "lightskyblue1", "lightskyblue2", "lightskyblue3", "lightskyblue4")))
export_figures(Fig16_The_Heat_Map_sign)


## Box plots on significant proteins:

Fig17_BoxPlots_sign <- GCP_BoxPlots(GCPlist = GCPlist14f,
                                    by_samples = FALSE)
export_figures(Fig17_BoxPlots_sign)                            ## exporting each box-plot as a single png file
export_figures(Fig17_BoxPlots_sign, exprt_fig_type = "pdf")    ## exporting all box-plots in a pdf file with many pages





## 4 mock groups for an ANOVA and Multiple FC analyes:

GCPlist14a <- GCPlist14
GCPlist14a$sampleINFO[, "Group_multi"] <- factor(c("gr1", "gr1", "gr2", "gr2", "gr3", "gr3", "gr3", "gr4", "gr4", "gr4"), levels = c("gr1", "gr2", "gr3", "gr4"))


GCPlist14a1 <- GCP_ANOVA(GCPlist = GCPlist14a,
                         name_column_groups = "Group_multi")


GCPlist14a2 <- GCP_FoldChange_Multi(GCPlist = GCPlist14a1,
                                    name_column_groups = "Group_multi",
                                    paired = FALSE,
                                    are_log_transf = TRUE)



# Batch correction!

GCPlist14b <- GCPlist14

## passing a numeric/factor vector directly to the 'batch' argument:

GCPlist14b1 <- GCP_ComBat(GCPlist = GCPlist14b,
                          batch = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1))

## passing a column name of sampleINFO, containing the numeric/factor vector, to the 'batch' argument:

GCPlist14b$sampleINFO[, "batches"] <- factor(c("batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1", "batch2", "batch3", "batch1"), levels = c("batch1", "batch2", "batch3"))

GCPlist14b2 <- GCP_ComBat(GCPlist = GCPlist14b,
                          batch = "batches")


# Exporting

## exporting a single GCPlist:

ExportGCPlist(GCPlist14, "GCPlist14_single_export.txt")


## exporting multiple GCPlists:

GCPlist_of_lists <- list(asimported = GCPlist00,
                         `before processing` = GCPlist05,
                         `post processing` = GCPlist11,
                         with_statistics = GCPlist14,
                         with_mock_ANOVA = GCPlist14a2,
                         with_mock_batchcorr = GCPlist14b2)

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
                                    GCPlistProteins = GCPlist14,
                                    are_log_transf = TRUE)

GCPlist_of_PTMlists <- list(asimported = GCPlistPTM00,
                            `before processing` = GCPlistPTM04,
                            `post processing` = GCPlistPTM08,
                            with_ratio_protein = GCPlistPTM09)
ExportGCPlist(GCPlist_of_PTMlists, "GCPlistsPTM_multiple_export.xlsx")

