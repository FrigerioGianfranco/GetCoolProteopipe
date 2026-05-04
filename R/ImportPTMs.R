#' Importing PTMs
#'
#' Starting from the output of PTMs, it stores the data in a list as the GCPlist.
#'
#' @param MaxQuant_table_name a character vector of length 1 with the name of the MaxQuant table file in the current working directory, which must be in the .txt format.
#' @param samples_info NULL or NA or a character vector of length 1 with the name of the table in the current working directory, containing information for each sample. The table must be in txt, csv, or xslsx format. In particular, the first column of the table must contain the names of the samples exactly as they are in the MaxQuant table.
#' @param fasta_database NULL or NA or either "human" or "mouse", or a name of a table in the current working directory (in .txt or .csv format). You can specify the fasta database to use to fill the missing Protein names and Gene names from the MaxQuant table. If you choose "mouse" or "human", the fasta table implemented were downloaded and reprocessed from Mascot on 5 May 2024.
#' @param prioritize_MaxQuant_names logical. If TRUE and if a fasta_database is provided, the final "Protein names" and "Gene names" will be primarily taken from the  MaxQuant table (they will be taken from the fasta database only if missing). If FALSE, the opposite will happen.
#' @param pattern_intensity a character vector of length 1. The character pattern that must be uniquely present in the column names of protein intensities.
#' @param pattern_isoforms NULL or a character vector of length 1. If not NULL, this function will import the different isoforms of each peptide as separate rows identified with the specified pattern (as example: "___1", "___2", "___3").
#'
#' @return a GCPlist, i.e.: a list with 3 dataframes (tibbles):
#'
#' - `intensities`: the intensities.
#'
#' - `proteinINFO`: all the information for each rows.
#'
#' - `sampleINFO`: the information about the samples.
#'
#'
#'
#' @examples
#' \dontrun{
#'
#'
#' GCPlistPTM00 <- ImportPTMs(MaxQuant_table_name = "PTM FILE NAME.txt",      ## you can just put the name of your file in your current working directory
#'                            samples_info = NULL,                             ## also here potentially
#'                            fasta_database = "mouse",
#'                            prioritize_MaxQuant_names = TRUE,
#'                            pattern_intensity = "Intensity ",
#'                            pattern_isoforms = "___")
#'
#' }
#'
#'
#' @importFrom readxl read_excel
#'
#' @export
ImportPTMs <- function(MaxQuant_table_name, samples_info = NULL, fasta_database = NA, prioritize_MaxQuant_names = TRUE, pattern_intensity = "Intensity ", pattern_isoforms = "___") {

  if (length(MaxQuant_table_name) != 1) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}
  if (is.na(MaxQuant_table_name)) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}
  if (!is.character(MaxQuant_table_name)) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}

  if (!is.null(samples_info)) {
    if (length(samples_info) != 1) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}

    if (!is.na(samples_info)) {
      if (!is.character(samples_info)) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}
    }
  }

  if (!is.null(fasta_database)) {

    if (length(fasta_database)!=1) {stop("if not NULL or NA, fasta_database must be a character of length 1")}

    if (!is.na(fasta_database)) {

      if (!is.character(fasta_database)) {stop("fasta_database must be a character")}

      if (tolower(fasta_database) == "na") {fasta_database <- NA}

      if (!is.na(fasta_database)) {

        if (tolower(fasta_database) == "human") {fasta_database <- "human"}
        if (tolower(fasta_database) == "mouse") {fasta_database <- "mouse"}

        if (fasta_database != "human" & fasta_database != "mouse" & !endsWith(fasta_database, ".txt") & !endsWith(fasta_database, ".csv")) {stop('if not NULL or NA, fasta_database must be "human", "mouse", or a name of a .txt or .csv table in the current working directory')}
      }
    }
  }


  if (!is.logical(prioritize_MaxQuant_names)) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}
  if (length(prioritize_MaxQuant_names) != 1) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}
  if (is.na(prioritize_MaxQuant_names)) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}

  if (length(pattern_intensity) != 1) {stop("pattern_intensity must be a character vector of length 1")}
  if (is.na(pattern_intensity)) {stop("pattern_intensity must be a character vector")}
  if (!is.character(pattern_intensity)) {stop("pattern_intensity must be a character vector")}

  if (!is.null(pattern_isoforms)) {
    if (length(pattern_isoforms) != 1) {stop("pattern_isoforms must be NULL or a character vector of length 1")}
    if (is.na(pattern_isoforms)) {stop("pattern_isoforms must be NULL or a character vector")}
    if (!is.character(pattern_isoforms)) {stop("pattern_isoforms must be NULL or a character vector")}
  }


  proteinGroup_Raw <- read_tsv(MaxQuant_table_name, guess_max = Inf)

  if(!(all(c("Protein", "Protein names", "Gene names", "id") %in% colnames(proteinGroup_Raw)))) {stop('MaxQuant_table must have at least the columns named: "Protein", "Protein names", "Gene names", "id"')}

  if (!is.character(proteinGroup_Raw$`Protein`)) stop('The column "Protein" must contain character')

  if (!is.character(proteinGroup_Raw$`Protein names`)) stop('The column "Protein names" must contain character')

  if (!is.character(proteinGroup_Raw$`Gene names`)) stop('The column "Gene names" must contain character')


  if (!is.numeric(proteinGroup_Raw$id)) stop('The column "id" must contain numbers')

  column_intesities_names <- colnames(proteinGroup_Raw)[grepl(pattern_intensity, colnames(proteinGroup_Raw))]

  if (length(column_intesities_names)<1) {stop("no coulumns with the specified pattern_intensity!")}


  if (is.null(pattern_isoforms)) {
    cat("\n___\n")
    cat(paste0("Considering the pattern_intensity '", pattern_intensity, "', the coulmns considered for the intensiteis are:\n '"))
    cat(paste0(column_intesities_names, collapse = "'\n '"))
    cat("'\n___\n")
  } else {
    if (any(grepl(pattern_isoforms, column_intesities_names))) {

      if ("isoform_indication"%in%colnames(proteinGroup_Raw)) {
        warning("there was already a column named 'isoform_indication'. Please, note that now it has been replaced")
        proteinGroup_Raw <- proteinGroup_Raw[, colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)!="isoform_indication")]]
      }

      table_separation_isoform <- tibble(original_name = column_intesities_names[which(grepl("___", column_intesities_names))],
                                         left_part = sub("___.*$", "", column_intesities_names[which(grepl("___", column_intesities_names))]),
                                         rigth_part = sub("^.*(?=___)", "", column_intesities_names[which(grepl("___", column_intesities_names))], perl = TRUE))
      if (any(duplicated(table_separation_isoform$original_name))) {stop("there are duplicated in samples name considering the chosen pattern_intensity and pattern_isoforms")}

      table_intensities_to_use <- proteinGroup_Raw[, table_separation_isoform$original_name]

      other_columns <- colnames(proteinGroup_Raw)[which(!colnames(proteinGroup_Raw)%in%table_separation_isoform$original_name)]

      new_intensities_col_names <- unique(table_separation_isoform$left_part)

      if (any(other_columns %in% new_intensities_col_names)) {

        col_name_intensity_tot <- other_columns[which(other_columns %in% new_intensities_col_names)]
        col_name_intensity_tot_with_tot <- paste0(col_name_intensity_tot, "_tot")

        if (any(col_name_intensity_tot_with_tot%in%colnames(proteinGroup_Raw))) {
          warning("there were column names already ending in '_tot', please not that they are being replaced to create columns that indicate the total intensities of a set of isoforms")
          proteinGroup_Raw <- proteinGroup_Raw[, colnames(proteinGroup_Raw)[which(!colnames(proteinGroup_Raw)%in%col_name_intensity_tot_with_tot)]]
        }

        colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)%in%col_name_intensity_tot)] <- paste0(colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)%in%col_name_intensity_tot)], "_tot")

        other_columns <- colnames(proteinGroup_Raw)[which(!colnames(proteinGroup_Raw)%in%table_separation_isoform$original_name)]
      }



      proteinGroup_Raw <- proteinGroup_Raw[, other_columns]
      proteinGroup_Raw <- mutate(proteinGroup_Raw, isoform_indication = rep(NA_character_, nrow(proteinGroup_Raw)))
      proteinGroup_Raw[, new_intensities_col_names] <- NA_character_
      proteinGroup_Raw_to_repeat <- proteinGroup_Raw


      for (single_isoform in unique(table_separation_isoform$rigth_part)) {

        table_separation_isoform_fil <- filter(table_separation_isoform, rigth_part == single_isoform)

        table_intensities_to_use_now <- table_intensities_to_use[, colnames(table_intensities_to_use)[which(colnames(table_intensities_to_use)%in%table_separation_isoform_fil$original_name)]]

        for (icaun in seq(colnames(table_intensities_to_use_now))) {
          colnames(table_intensities_to_use_now)[icaun] <- table_separation_isoform_fil$left_part[which(table_separation_isoform_fil$original_name==colnames(table_intensities_to_use_now)[icaun])]
        }


        if (single_isoform == unique(table_separation_isoform$rigth_part)[1]) {

          proteinGroup_Raw$isoform_indication <- rep(single_isoform, nrow(proteinGroup_Raw))
          for (cntun in colnames(table_intensities_to_use_now)) {
            proteinGroup_Raw[, cntun] <- pull(table_intensities_to_use_now, cntun)
          }

        } else {

          df_to_bindrow_below <- proteinGroup_Raw_to_repeat
          df_to_bindrow_below$isoform_indication <- rep(single_isoform, nrow(df_to_bindrow_below))
          for (cntun in colnames(table_intensities_to_use_now)) {
            df_to_bindrow_below[, cntun] <- pull(table_intensities_to_use_now, cntun)
          }
          proteinGroup_Raw <- bind_rows(proteinGroup_Raw, df_to_bindrow_below)
        }
      }

      cat("\n___\n")
      cat(paste0("Considering the pattern_intensity '", pattern_intensity, "', and the pattern_isoforms '", pattern_isoforms, "' the following sample intensities were considered:\n '"))
      cat(paste0(new_intensities_col_names, collapse = "'\n '"))
      cat("'\n\nrepeated for the following isoforms:\n '")
      cat(paste0(unique(table_separation_isoform$rigth_part), collapse = "'\n '"))
      cat("'\n\n___\n")


      column_intesities_names <- new_intensities_col_names

    } else {
      stop(paste0("The pattern_isoforms '", pattern_isoforms, "' was not found in any of the column names idendified with the pattern_intensity '", pattern_intensity, "'"))
    }
  }


  for (a in column_intesities_names) {
    if (!is.numeric(pull(proteinGroup_Raw, a))) {
      proteinGroup_Raw[,a] <- as.numeric(pull(proteinGroup_Raw, a))
    }
  }

  if ("protid" %in% colnames(proteinGroup_Raw)) {
    warning("the MaxQuant table already had a column named 'protid', now it has been replaced")
    proteinGroup_Raw <- proteinGroup_Raw[, colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)!="protid")]]
  }


  samplenames_raw_intensities <- column_intesities_names

  samplenames_raw_clean <- samplenames_raw_intensities

  if (any(duplicated(samplenames_raw_clean))) {stop("The names of samples must not be duplicated. Each of them must be unique!")}


  if (!is.null(samples_info)) {
    if (!is.na(samples_info)) {
      if (endsWith(toupper(samples_info), ".CSV")) {
        samples_info_loaded <- read_csv(samples_info)
      } else if (endsWith(toupper(samples_info), ".XLSX")) {
        samples_info_loaded <- read_excel(samples_info, sheet = 1)
      } else {
        samples_info_loaded <- read_tsv(samples_info)
      }


      if (!all(pull(samples_info_loaded, 1) %in% samplenames_raw_clean)) {
        stop("The names contained in the first column of the samples_info table are not matching with the sample names in the MaxQuant table")
      }

      if (any(duplicated(pull(samples_info_loaded, 1)))) {stop("The names of samples from the first column of samples_info must not contain duplicate. Each of them must be unique")}

      for (a in colnames(samples_info_loaded)) {
        if (is.character(pull(samples_info_loaded, a))) {
          if (any(duplicated(pull(samples_info_loaded, a)))) {
            samples_info_loaded[,a] <- as.factor(pull(samples_info_loaded, a))
          }
        }
      }


    } else {
      samples_info_loaded <- tibble(Sample = samplenames_raw_clean)
    }
  } else {
    samples_info_loaded <- tibble(Sample = samplenames_raw_clean)
  }




  if (!is.null(fasta_database)) {
    if (!is.na(fasta_database)) {
      if (fasta_database == "human") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Human_ref_20240305.txt.gz", package = "GetCoolProteopipe"), show_col_types = FALSE)
      } else if (fasta_database == "mouse") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Mouse_ref_20240305.txt.gz", package = "GetCoolProteopipe"), show_col_types = FALSE)
      } else if (endsWith(fasta_database, ".txt")) {
        fasta_database_loaded <- read_tsv(fasta_database, show_col_types = FALSE)
      } else if (endsWith(fasta_database, ".csv")) {
        fasta_database_loaded <- read_csv(fasta_database, show_col_types = FALSE)
      }

      if (!(all(c("Accession", "Protein names", "Gene names") %in% colnames(fasta_database_loaded))))  {stop('The fasta database loaded must contain at least these columns: "Accession", "Protein names", and "Gene names"')}
      if (any(duplicated(pull(fasta_database_loaded, "Accession")))) {stop('The fasta database must not contain any duplicated in the "Accession" column')}
    }
  }


  if (any(duplicated(proteinGroup_Raw$id))) {

    if (!is.null(pattern_isoforms)) {
      if (all(startsWith(proteinGroup_Raw$isoform_indication, "_"))) {
        id_deduplicated <- paste0(as.character(proteinGroup_Raw$id), proteinGroup_Raw$isoform_indication)
      } else {
        id_deduplicated <- paste0(as.character(proteinGroup_Raw$id), rep("_", nrow(proteinGroup_Raw)), proteinGroup_Raw$isoform_indication)
      }

      if (any(duplicated(id_deduplicated))) {
        id_deduplicated <- fix_duplicated(id_deduplicated,
                                          zeros = TRUE,
                                          define_highest_for_zeros = NULL,
                                          start_with_zero = FALSE,
                                          exclude_the_first = FALSE,
                                          NA_as_character = TRUE)
      }
    } else {
      id_deduplicated <- fix_duplicated(as.character(proteinGroup_Raw$id),
                                        zeros = TRUE,
                                        define_highest_for_zeros = NULL,
                                        start_with_zero = FALSE,
                                        exclude_the_first = FALSE,
                                        NA_as_character = TRUE)
    }
  } else {
    id_deduplicated <- as.character(proteinGroup_Raw$id)
  }


  proteinGroup_Raw <- proteinGroup_Raw %>%
    add_column(protid = "PTM",
               .before = 1) %>%
    mutate(protid = ifelse(id<10, paste0("PTM00000", id_deduplicated),
                         ifelse(id<100, paste0("PTM0000", id_deduplicated),
                                ifelse(id<1000, paste0("PTM000", id_deduplicated),
                                       ifelse(id<10000, paste0("PTM00",id_deduplicated),
                                              ifelse(id<100000, paste0("PTM0", id_deduplicated),
                                                     paste0("PTM", id_deduplicated)))))))



  if ("Accession" %in% colnames(proteinGroup_Raw)) {warning('A column named "Accession" was already present in the imported table. Please not that now it has been completely replaced with each first code (before the ";") of "Protein"')}

  proteinGroup_Raw_acc <-  mutate(proteinGroup_Raw, Accession = gsub(";.*", "", `Protein`))

  if (!is.null(fasta_database)) {
    if (!is.na(fasta_database)) {

      proteinGroup_Raw_fastajoined <- left_join(x = proteinGroup_Raw_acc, y = fasta_database_loaded, by = "Accession", suffix = c("_MaxQuant", "_fasta"))

      if (prioritize_MaxQuant_names) {
        proteinGroup_Raw_added <- mutate(proteinGroup_Raw_fastajoined,
                                         `Protein names` = ifelse(is.na(`Protein names_MaxQuant`), `Protein names_fasta`, `Protein names_MaxQuant`),
                                         `Gene names` = ifelse(is.na(`Gene names_MaxQuant`), `Gene names_fasta`, `Gene names_MaxQuant`))

      } else {

        proteinGroup_Raw_added <- mutate(proteinGroup_Raw_fastajoined,
                                         `Protein names` = ifelse(is.na(`Protein names_fasta`), `Protein names_MaxQuant`, `Protein names_fasta`),
                                         `Gene names` = ifelse(is.na(`Gene names_fasta`), `Gene names_MaxQuant`, `Gene names_fasta`))


      }
    } else {
      proteinGroup_Raw_added <- proteinGroup_Raw_acc
      colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Protein names")] <- "Protein names_MaxQuant"
      colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Gene names")] <- "Gene names_MaxQuant"

      proteinGroup_Raw_added <- mutate(proteinGroup_Raw_added,
                                       `Protein names` = `Protein names_MaxQuant`,
                                       `Gene names` = `Gene names_MaxQuant`)
    }
  } else {
    proteinGroup_Raw_added <- proteinGroup_Raw_acc
    colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Protein names")] <- "Protein names_MaxQuant"
    colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Gene names")] <- "Gene names_MaxQuant"

    proteinGroup_Raw_added <- mutate(proteinGroup_Raw_added,
                                     `Protein names` = `Protein names_MaxQuant`,
                                     `Gene names` = `Gene names_MaxQuant`)
  }


  proteinGroup_Raw_added_cut <- mutate(proteinGroup_Raw_added,
                                       `Protein names` = gsub(";.*", "", `Protein names`),
                                       `Gene names` = gsub(";.*", "", `Gene names`))


  colnames_raw_intensities <- c("protid", samplenames_raw_intensities)

  colnames_others <- c("protid", colnames(proteinGroup_Raw_added_cut)[which(!(colnames(proteinGroup_Raw_added_cut)%in%colnames_raw_intensities[-1]) & !(colnames(proteinGroup_Raw_added_cut)=="protid"))])


  final_list <- list(intensities = select(proteinGroup_Raw_added_cut, all_of(colnames_raw_intensities)),
                     proteinINFO = select(proteinGroup_Raw_added_cut, all_of(colnames_others)),
                     sampleINFO = samples_info_loaded)

  final_list$intensities <- mutate_at(final_list$intensities,
                                    colnames(final_list$intensities)[which(colnames(final_list$intensities) != "protid")],
                                    ~ replace(., is.nan(.), NA_real_))

  return(final_list)

}

