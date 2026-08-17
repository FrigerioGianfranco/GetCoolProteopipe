#' Importing the output of Proteome Discoverer
#'
#' Starting from the output of the Proteome Discoverer, it performs some cleaning such as filtering out proteins that are only identified by site, the reverse and the potential contaminants, then it adds protein names from a fasta database.
#'
#' @param ProtDiscov_table_name a character vector of length 1 with the name of the Proteome Discoverer table file, as exported, in the current working directory, which must be in the .txt format.
#' @param ProtDiscov_InputFiles  NULL or a character vector of length 1 with the name of the Proteome Discoverer input file table, as exported, in the current working directory, which must be in the .txt format.
#' @param samples_info NULL or NA or a character vector of length 1 with the name of the table in the current working directory, containing information for each sample. The table must be in txt, csv, or xslsx format. In particular, the first column of the table must contain the names of the samples exactly the same considered.
#' @param restore_sample_names logical. If TRUE, you must provide the ProtDiscov_InputFiles and the original sample names will be considered. If FALSE, the Proteome Discoverer sample names (like F1, F2, ...) will be considered.
#' @param remove_empty_columns logical. If TRUE, samples that have full missing values abundances will be removed.
#' @param fasta_database NULL or NA or either "human" or "mouse", or a name of a table in the current working directory (in .txt or .csv format). You can specify the fasta database to use to fill the missing Protein names and Gene names from the Proteome Discoverer table. If you choose "mouse" or "human", the fasta table implemented were downloaded and reprocessed from Mascot on 5 May 2024.
#' @param prioritize_ProtDiscov_names logical. If TRUE and if a fasta_database is provided, the final "Protein names" and "Gene names" will be primarily taken from the  Proteome Discoverer table (they will be taken from the fasta database only if missing). If FALSE, the opposite will happen.
#' @param remove_contaminants logical. Do you want to remove rows that are TRUE in the column "Contaminant"?
#'
#' @return a GCPlist, i.e.: a list with 3 dataframes (tibbles):
#'
#' - `intensities`: the protein intensities.
#'
#' - `proteinINFO`: all the information for each rows.
#'
#' - `sampleINFO`: the information about the samples.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist00P <- ImportOutputProtDiscov(ProtDiscov_table_name = "Proteome Discoverer OUTPUT FILE NAME.txt",     ## put here the name of the Proteome Discoverer output table present in your current working directory
#'                                      ProtDiscov_InputFiles = "Proteome Discoverer INPUT FILE NAME.txt",      ## you could also put the name the input file if you want to restore the original sample names
#'                                      samples_info = NULL,                                                 ## you could put here the file name of a table with further information about your samples
#'                                      restore_sample_names = TRUE,
#'                                      remove_empty_columns = TRUE,
#'                                      fasta_database = "human",
#'                                      prioritize_ProtDiscov_names = TRUE,
#'                                      remove_contaminants = TRUE)
#'
#' }
#'
#' @importFrom readxl read_excel
#'
#' @export
ImportOutputProtDiscov <- function(ProtDiscov_table_name, ProtDiscov_InputFiles = NULL, samples_info = NULL, restore_sample_names = TRUE, remove_empty_columns = TRUE, fasta_database = NA, prioritize_ProtDiscov_names = TRUE, remove_contaminants = TRUE) {

  if (length(ProtDiscov_table_name) != 1) {stop("ProtDiscov_table_name must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}
  if (is.na(ProtDiscov_table_name)) {stop("ProtDiscov_table_name must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}
  if (!is.character(ProtDiscov_table_name)) {stop("ProtDiscov_table_name must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}

  if (!is.null(ProtDiscov_InputFiles)) {
    if (length(ProtDiscov_InputFiles) != 1) {stop("ProtDiscov_InputFiles must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}
    if (is.na(ProtDiscov_InputFiles)) {stop("ProtDiscov_InputFiles must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}
    if (!is.character(ProtDiscov_InputFiles)) {stop("ProtDiscov_InputFiles must be a character vector of length 1, indicating the name of the Proteome Discoverer table files, in txt format")}
  }

  if (!is.null(samples_info)) {
    if (length(samples_info) != 1) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}

    if (!is.na(samples_info)) {
      if (!is.character(samples_info)) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}
    }
  }

  if (!is.logical(restore_sample_names)) {stop("restore_sample_names must be either TRUE or FALSE")}
  if (length(restore_sample_names) != 1) {stop("restore_sample_names must be either TRUE or FALSE")}
  if (is.na(restore_sample_names)) {stop("restore_sample_names must be either TRUE or FALSE")}
  if (restore_sample_names & is.null(ProtDiscov_InputFiles)) {stop("Since restore_sample_names is TRUE, you must provide the ProtDiscov_InputFiles. Otherwise, set restore_sample_names as FALSE")}

  if (!is.logical(remove_empty_columns)) {stop("remove_empty_columns must be either TRUE or FALSE")}
  if (length(remove_empty_columns) != 1) {stop("remove_empty_columns must be either TRUE or FALSE")}
  if (is.na(remove_empty_columns)) {stop("remove_empty_columns must be either TRUE or FALSE")}

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

  if (!is.logical(prioritize_ProtDiscov_names)) {stop("prioritize_ProtDiscov_names must be either TRUE or FALSE")}
  if (length(prioritize_ProtDiscov_names) != 1) {stop("prioritize_ProtDiscov_names must be either TRUE or FALSE")}
  if (is.na(prioritize_ProtDiscov_names)) {stop("prioritize_ProtDiscov_names must be either TRUE or FALSE")}

  if (!is.logical(remove_contaminants)) {stop("remove_contaminants must be either TRUE or FALSE")}
  if (length(remove_contaminants) != 1) {stop("remove_contaminants must be either TRUE or FALSE")}
  if (is.na(remove_contaminants)) {stop("remove_contaminants must be either TRUE or FALSE")}


  proteinGroup_Raw <- read_tsv(ProtDiscov_table_name, guess_max = Inf)
  if (!is.null(ProtDiscov_InputFiles)) {
    PD_InputFiles_tab <- read_tsv(ProtDiscov_InputFiles, guess_max = Inf)
    if (!"File Name"%in%colnames(PD_InputFiles_tab)) {stop("'File Name' must be present as a column in the ProtDiscov_InputFiles table")}
    if (!"File ID"%in%colnames(PD_InputFiles_tab)) {stop("'File ID' must be present as a column in the ProtDiscov_InputFiles table")}
    if (any(duplicated(PD_InputFiles_tab$`File ID`))) {stop("There are duplicated in the 'File ID' column of the ProtDiscov_InputFiles table")}
  }

  if(!(all(c("Accession", "Description", "Gene Symbol") %in% colnames(proteinGroup_Raw)))) {stop('ProtDiscov_table must have at least the columns named: "Accession", "Description", "Gene Symbol"')}
  if (!is.character(proteinGroup_Raw$`Accession`)) {stop('The column "Accession" must contain character')}
  if (!is.character(proteinGroup_Raw$`Description`)) {stop('The column "Description" must contain character')}
  if (!is.character(proteinGroup_Raw$`Gene Symbol`)) {stop('The column "Gene Symbol" must contain character')}

  if (remove_contaminants) {
    if (!("Contaminant" %in% colnames(proteinGroup_Raw))) {stop("'Contaminant' must be a column of the ProtDiscov_table\n\n ...you could also set the argument remove_contaminants as FALSE to avoid this error!")}
    if (!is.logical(proteinGroup_Raw$Contaminant)) {stop("'Contaminant' in the ProtDiscov_table must contain only TRUE/FALSE")}
    if (any(is.na(proteinGroup_Raw$Contaminant))) {stop("'Contaminant' in the ProtDiscov_table must not contain missing values")}
  }

  if ("Protein names" %in% colnames(proteinGroup_Raw)) {
    warning("A column named 'Protein names' was already present in the ProtDiscov table, but now it has been replaced with the content of the 'Description' column")
    proteinGroup_Raw <- proteinGroup_Raw[, which(colnames(proteinGroup_Raw)!="Protein names")]
  }
  colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)=="Description")] <- "Protein names"

  if ("Gene names" %in% colnames(proteinGroup_Raw)) {
    warning("A column named 'Gene names' was already present in the ProtDiscov table, but now it has been replaced with the content of the 'Gene Symbol' column")
    proteinGroup_Raw <- proteinGroup_Raw[, which(colnames(proteinGroup_Raw)!="Gene names")]
  }
  colnames(proteinGroup_Raw)[which(colnames(proteinGroup_Raw)=="Gene Symbol")] <- "Gene names"


  colnames_with_intensities <- colnames(proteinGroup_Raw)[which(grepl("Abundance", colnames(proteinGroup_Raw)))]
  colnames_with_averaged_intensities <- colnames(proteinGroup_Raw)[which(grepl("Abundances Grouped", colnames(proteinGroup_Raw)))]
  colnames_with_sample_intensities <- colnames_with_intensities[which(!colnames_with_intensities%in%colnames_with_averaged_intensities)]

  if (length(colnames_with_intensities) < 1) {stop("There are no 'Abundance' columns in the Proteome Discoverer table!")}
  if (length(colnames_with_sample_intensities) < 1) {stop("There are no 'Abundance' columns related to specific samples in the Proteome Discoverer table besides the 'Abundances Grouped'!")}

  sample_intensities_names_df <- tibble(original_colnames_with_sample_intensities = colnames_with_sample_intensities,
                                        names_from_PD = sub("^Abundance\\s+(.*?)\\s+Sample.*$", "\\1", colnames_with_sample_intensities),
                                        groups_from_PD = sub("^.*Sample\\s+", "", colnames_with_sample_intensities),
                                        has_all_NA = rep(FALSE, length(colnames_with_sample_intensities)))
  if (any(duplicated(sample_intensities_names_df$original_colnames_with_sample_intensities))) {stop("There are duplicated in the Abundance columns of the Proteome Discoverer table..")}
  if (any(duplicated(sample_intensities_names_df$names_from_PD))) {stop("There are duplicated in the Abundance column sample names of the Proteome Discoverer table..")}

  real_sample_names_used <- FALSE

  if (!is.null(ProtDiscov_InputFiles)) {

    file_pattern_to_find <- "^.*\\\\([^\\\\]+)\\.raw$"

    names_from_File_Name_vctr <- ifelse(grepl(file_pattern_to_find, pull(PD_InputFiles_tab, "File Name")),
                                        sub(file_pattern_to_find, "\\1", pull(PD_InputFiles_tab, "File Name")),
                                        NA_character_)

    sample_intensities_names_df <- add_column(sample_intensities_names_df,
                                              names_from_File_Name = rep(NA_character_, nrow(sample_intensities_names_df)))
    for (i_sindf in seq(nrow(sample_intensities_names_df))) {
      if (sample_intensities_names_df$names_from_PD[i_sindf] %in% PD_InputFiles_tab$`File ID`) {
        sample_intensities_names_df$names_from_File_Name[i_sindf] <- names_from_File_Name_vctr[which(PD_InputFiles_tab$`File ID`==sample_intensities_names_df$names_from_PD[i_sindf])]
      }
    }

    if (restore_sample_names) {
      real_sample_names_used <- TRUE

      if (any(is.na(sample_intensities_names_df$names_from_File_Name))) {
        NA_info_to_print <- character()
        for (i_NA_orig_names in which(is.na(sample_intensities_names_df$names_from_File_Name))) {
          NA_info_to_print <- c(NA_info_to_print, paste0(sample_intensities_names_df$names_from_PD[i_NA_orig_names], " : ", sample_intensities_names_df$names_from_File_Name[i_NA_orig_names]))
        }
        stop(paste0("Sample names cannot be restored as there are names that cannot be matched or are NAs in the sample names from the InputFiles table:\n",
                    paste0(NA_info_to_print, collapse = "\n"),
                    "\n\n...You could also set the argument restore_sample_names as FALSE to avoid this error!"))
      }

      if (any(sample_intensities_names_df$names_from_File_Name == "")) {
        empty_info_to_print <- character()
        for (i_empty_orig_names in which(sample_intensities_names_df$names_from_File_Name == "")) {
          empty_info_to_print <- c(empty_info_to_print, paste0(sample_intensities_names_df$names_from_PD[i_empty_orig_names], " : ", sample_intensities_names_df$names_from_File_Name[i_empty_orig_names]))
        }
        stop(paste0("Sample names cannot be restored as there are empty names generated from the InputFiles table:\n",
                    paste0(empty_info_to_print, collapse = "\n"),
                    "\n\n...You could also set the argument restore_sample_names as FALSE to avoid this error!"))
      }

      if (any(duplicated(sample_intensities_names_df$names_from_File_Name))) {
        sample_intensities_names_df_dupl_index <- which(sample_intensities_names_df$names_from_File_Name%in%unique(sample_intensities_names_df$names_from_File_Name[which(duplicated(sample_intensities_names_df$names_from_File_Name))]))
        duplicated_info_to_print <- character()
        for (i_dupl in sample_intensities_names_df_dupl_index) {
          duplicated_info_to_print <- c(duplicated_info_to_print, paste0(sample_intensities_names_df$names_from_PD[i_dupl], " : ", sample_intensities_names_df$names_from_File_Name[i_dupl]))
        }
        stop(paste0("Sample names cannot be restored as there are duplicated in the sample names from the InputFiles table:\n",
                    paste0(duplicated_info_to_print, collapse = "\n"),
                    "\n\n...You could also set the argument restore_sample_names as FALSE to avoid this error!"))
      }
    }
  }

  if (any(map_lgl(select(proteinGroup_Raw, all_of(colnames_with_sample_intensities)), \(z) all(is.na(z))))) {
    for (iNA in which(colnames_with_sample_intensities %in% names(map_lgl(select(proteinGroup_Raw, all_of(colnames_with_sample_intensities)), \(z) all(is.na(z)))[which(map_lgl(select(proteinGroup_Raw, all_of(colnames_with_sample_intensities)), \(z) all(is.na(z))))]))) {
      sample_intensities_names_df[iNA, "has_all_NA"] <- TRUE
      proteinGroup_Raw[,colnames_with_sample_intensities[iNA]] <- as.numeric(pull(proteinGroup_Raw, colnames_with_sample_intensities[iNA]))
    }
  }

  if (!(all(map_lgl(select(proteinGroup_Raw, all_of(colnames_with_sample_intensities)), is.numeric)))) {
    stop("all the Abundance columns must contain numerical data")
  }

  if ("protid" %in% colnames(proteinGroup_Raw)) {
    warning("The Proteome Discoverer table contained already a column named 'protid', but now it has been completely replaced!")
    proteinGroup_Raw <- proteinGroup_Raw[, which(colnames(proteinGroup_Raw)!="protid")]
  }
  protid_vctr <- character()
  for (i_id in seq(from = 0, to = (nrow(proteinGroup_Raw)-1))) {
    protid_vctr <- c(protid_vctr, ifelse(i_id<10, paste0("prot00000", i_id),
                                         ifelse(i_id<100, paste0("prot0000", i_id),
                                                ifelse(i_id<1000, paste0("prot000", i_id),
                                                       ifelse(i_id<10000, paste0("prot00", i_id),
                                                              ifelse(i_id<100000, paste0("prot0", i_id),
                                                                     paste0("prot", i_id)))))))
  }
  proteinGroup_Raw <- add_column(proteinGroup_Raw,
                                 protid = protid_vctr,
                                 .before = 1)


  proteinGroup_Raw_acc <- proteinGroup_Raw[, which(!colnames(proteinGroup_Raw)%in%colnames_with_sample_intensities)]

  the_intensity_df <- proteinGroup_Raw[, c(1, which(colnames(proteinGroup_Raw)%in%colnames_with_sample_intensities))]

  if (!real_sample_names_used) {
    sample_names_used <- sample_intensities_names_df$names_from_PD

    for (id_col_intensity_df in seq(from = 2, to = length(colnames(the_intensity_df)))) {
      colnames(the_intensity_df)[id_col_intensity_df] <- sample_intensities_names_df$names_from_PD[which(sample_intensities_names_df$original_colnames_with_sample_intensities==colnames(the_intensity_df)[id_col_intensity_df])]
    }

    starting_sampleINFO <- tibble(Sample = sample_names_used,
                                  Groups_from_PD = as.factor(sample_intensities_names_df$groups_from_PD))
    if (!is.null(ProtDiscov_InputFiles)) {
      starting_sampleINFO <- add_column(starting_sampleINFO,
                                        names_from_File_Name = sample_intensities_names_df$names_from_File_Name,
                                        .after = 1)
    }
  } else {
    sample_names_used <- sample_intensities_names_df$names_from_File_Name

    for (id_col_intensity_df in seq(from = 2, to = length(colnames(the_intensity_df)))) {
      colnames(the_intensity_df)[id_col_intensity_df] <- sample_intensities_names_df$names_from_File_Name[which(sample_intensities_names_df$original_colnames_with_sample_intensities==colnames(the_intensity_df)[id_col_intensity_df])]
    }

    starting_sampleINFO <- tibble(Sample = sample_names_used,
                                  names_from_PD = sample_intensities_names_df$names_from_PD,
                                  Groups_from_PD = as.factor(sample_intensities_names_df$groups_from_PD))
  }


  if (!is.null(samples_info)) {
    if (!is.na(samples_info)) {
      if (endsWith(toupper(samples_info), ".CSV")) {
        samples_info_loaded <- read_csv(samples_info)
      } else if (endsWith(toupper(samples_info), ".XLSX")) {
        samples_info_loaded <- read_excel(samples_info, sheet = 1)
      } else {
        samples_info_loaded <- read_tsv(samples_info)
      }


      if (!all(pull(samples_info_loaded, 1) %in% sample_names_used)) {
        stop("The names contained in the first column of the samples_info table are not matching with the sample names we are considering here\n\n...You could also leave the argument samples_info as NULL to avoid this error!")
      }

      if (any(duplicated(pull(samples_info_loaded, 1)))) {stop("The names of samples from the first column of samples_info must not contain duplicate")}

      for (a in colnames(samples_info_loaded)) {
        if (is.character(pull(samples_info_loaded, a))) {
          if (any(duplicated(pull(samples_info_loaded, a)))) {
            samples_info_loaded[,a] <- as.factor(pull(samples_info_loaded, a))
          }
        }
      }

      colnames(starting_sampleINFO)[1] <- colnames(samples_info_loaded)[1]
      if (any(duplicated(colnames(starting_sampleINFO)))) {
        colnames(starting_sampleINFO) <- fix_duplicated(colnames(starting_sampleINFO))
      }

      samples_info_considered <- left_join(x = starting_sampleINFO, y = samples_info_loaded, by = colnames(samples_info_loaded)[1], suffix = c("_PD", "_loaded"))

    } else {
      samples_info_considered <- starting_sampleINFO
    }
  } else {
    samples_info_considered <- starting_sampleINFO
  }


  if (!is.null(fasta_database)) {
    if (!is.na(fasta_database)) {
      if (fasta_database == "human") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Human_ref_20240305.txt.gz", package = "GetCoolProteopipe"))
      } else if (fasta_database == "mouse") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Mouse_ref_20240305.txt.gz", package = "GetCoolProteopipe"))
      } else if (endsWith(fasta_database, ".txt")) {
        fasta_database_loaded <- read_tsv(fasta_database)
      } else if (endsWith(fasta_database, ".csv")) {
        fasta_database_loaded <- read_csv(fasta_database)
      }

      if (!(all(c("Accession", "Protein names", "Gene names") %in% colnames(fasta_database_loaded))))  {stop('The fasta database loaded must contain at least these columns: "Accession", "Protein names", and "Gene names"')}
      if (any(duplicated(pull(fasta_database_loaded, "Accession")))) {stop('The fasta database must not contain any duplicated in the "Accession" column')}
    }
  }


  cat("\n___\n")

  if (real_sample_names_used) {
    cat(paste0('Original sample names successfully restored!\n (For example:\n  "',
               sample_intensities_names_df$names_from_PD[1], '" is "', sample_intensities_names_df$names_from_File_Name[1], '"\n'))
    if (nrow(sample_intensities_names_df)>=2) {cat(paste0('  "', sample_intensities_names_df$names_from_PD[2], '" is "', sample_intensities_names_df$names_from_File_Name[2], '"\n'))}
    if (nrow(sample_intensities_names_df)>=3) {cat(paste0('  "', sample_intensities_names_df$names_from_PD[3], '" is "', sample_intensities_names_df$names_from_File_Name[3], '"\n'))}
    cat('  and so on...)\n___\n')
  }

  if (remove_empty_columns) {
    if (any(sample_intensities_names_df$has_all_NA)) {

      cat("The following samples contained only missing values, so they have been removed:\n")

      if (real_sample_names_used) {
        samples_about_to_be_removed <- sample_intensities_names_df$names_from_File_Name[which(sample_intensities_names_df$has_all_NA)]

        samples_about_to_be_removed_to_print <- character()

        for (i_all_NA in which(sample_intensities_names_df$has_all_NA)) {
          samples_about_to_be_removed_to_print <- c(samples_about_to_be_removed_to_print, paste0(' "', sample_intensities_names_df$names_from_PD[i_all_NA], '": "', sample_intensities_names_df$names_from_File_Name[i_all_NA], '"'))
        }
      } else {
        samples_about_to_be_removed <- sample_intensities_names_df$names_from_PD[which(sample_intensities_names_df$has_all_NA)]

        samples_about_to_be_removed_to_print <- character()

        if (!is.null(ProtDiscov_InputFiles)) {
          for (i_all_NA in which(sample_intensities_names_df$has_all_NA)) {
            samples_about_to_be_removed_to_print <- c(samples_about_to_be_removed_to_print, paste0(' "', sample_intensities_names_df$names_from_PD[i_all_NA], '": "', sample_intensities_names_df$names_from_File_Name[i_all_NA], '"'))
          }
        } else {
          for (i_all_NA in which(sample_intensities_names_df$has_all_NA)) {
            samples_about_to_be_removed_to_print <- c(samples_about_to_be_removed_to_print, paste0(' "', sample_intensities_names_df$names_from_PD[i_all_NA], '"'))
          }
        }
      }

      cat(paste0(samples_about_to_be_removed_to_print, collapse = "\n"))
      cat('\n___\n')

      samples_info_considered <- samples_info_considered[which(!pull(samples_info_considered, 1)%in%samples_about_to_be_removed),]

      the_intensity_df <- the_intensity_df[, which(!colnames(the_intensity_df)%in%samples_about_to_be_removed)]
    }
  }

  if (remove_contaminants) {
    cat(paste0("At first, ", length(pull(proteinGroup_Raw_acc, 1))  , " rows were imported\n",
               "...from which, ", sum(proteinGroup_Raw_acc$Contaminant), " contaminants were removed  (those rows were TRUE in the 'Contaminant' column)\n"))

    proteinGroup_Raw_acc <- filter(proteinGroup_Raw_acc, !Contaminant)

    cat(paste0("In the end, ", length(pull(proteinGroup_Raw_acc, 1))  , " rows were kept"))

  } else {
    cat(paste0(length(pull(proteinGroup_Raw_acc, 1))  , " rows were imported and kept"))
  }
  cat("\n___\n")

  proteinGroup_Raw_acc <-  mutate(proteinGroup_Raw_acc, Accession = gsub(";.*", "", Accession))

  if (!is.null(fasta_database)) {
    if (!is.na(fasta_database)) {

      proteinGroup_Raw_fastajoined <- left_join(x = proteinGroup_Raw_acc, y = fasta_database_loaded, by = "Accession", suffix = c("_ProtDiscov", "_fasta"))

      if (prioritize_ProtDiscov_names) {
        proteinGroup_Raw_added <- mutate(proteinGroup_Raw_fastajoined,
                                         `Protein names` = ifelse(is.na(`Protein names_ProtDiscov`), `Protein names_fasta`, `Protein names_ProtDiscov`),
                                         `Gene names` = ifelse(is.na(`Gene names_ProtDiscov`), `Gene names_fasta`, `Gene names_ProtDiscov`))

      } else {

        proteinGroup_Raw_added <- mutate(proteinGroup_Raw_fastajoined,
                                         `Protein names` = ifelse(is.na(`Protein names_fasta`), `Protein names_ProtDiscov`, `Protein names_fasta`),
                                         `Gene names` = ifelse(is.na(`Gene names_fasta`), `Gene names_ProtDiscov`, `Gene names_fasta`))


      }
    } else {
      proteinGroup_Raw_added <- proteinGroup_Raw_acc
      colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Protein names")] <- "Protein names_ProtDiscov"
      colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Gene names")] <- "Gene names_ProtDiscov"

      proteinGroup_Raw_added <- mutate(proteinGroup_Raw_added,
                                       `Protein names` = `Protein names_ProtDiscov`,
                                       `Gene names` = `Gene names_ProtDiscov`)
    }
  } else {
    proteinGroup_Raw_added <- proteinGroup_Raw_acc
    colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Protein names")] <- "Protein names_ProtDiscov"
    colnames(proteinGroup_Raw_added)[which(colnames(proteinGroup_Raw_added) == "Gene names")] <- "Gene names_ProtDiscov"

    proteinGroup_Raw_added <- mutate(proteinGroup_Raw_added,
                                     `Protein names` = `Protein names_ProtDiscov`,
                                     `Gene names` = `Gene names_ProtDiscov`)
  }


  proteinGroup_Raw_added_cut <- mutate(proteinGroup_Raw_added,
                                       `Protein names` = gsub(";.*", "", `Protein names`),
                                       `Gene names` = gsub(";.*", "", `Gene names`))

  final_list <- list(intensities = the_intensity_df[which(the_intensity_df$protid %in% proteinGroup_Raw_added_cut$protid),],
                     proteinINFO = proteinGroup_Raw_added_cut,
                     sampleINFO = samples_info_considered)

  return(final_list)
}

