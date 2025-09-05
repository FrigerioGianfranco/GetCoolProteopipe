#' Importing PTMs
#'
#' Starting from the output of PTMs, it stores the raw data in a list as the GCPlist.
#'
#' @param MaxQuant_table_name a character vector of length 1 with the name of the MaxQuant table file in the current working directory, which must be in the .txt format.
#' @param samples_info NULL or NA or a character vector of length 1 with the name of the table in the current working directory, containing information for each sample. The table must be in txt, csv, or xslsx format. In particular, the first column of the table must contain the names of the samples exactly as they are in the MaxQuant table.
#' @param fasta_database NULL or NA or either "human" or "mouse", or a name of a table in the current working directory (in .txt or .csv format). You can specify the fasta database to use to fill the missing Protein names and Gene names from the MaxQuant table.
#' @param prioritize_MaxQuant_names logical. If TRUE and if a fasta_database is provided, the final "Protein names" and "Gene names" will be primarily taken from the  MaxQuant table (they will be taken from the fasta database only if missing). If FALSE, the opposite will happen.
#' @param pattern_intensity a character vector of length 1. The character pattern that must be uniquely present in the column names of protein intensities.
#'
#' @return a GCPlist, i.e.: a list with 4 dataframes (tibbles):#'
#'
#' - `quant_raw`: the raw intensities.
#'
#' - `quant_LFQ`: the LFQ intensities. (all NA fo this import)
#'
#' - `proteinINFO`: all the information for each rows.
#'
#' - `sampleINFO`: the information about the samples.
#'
#'
#' @importFrom readxl read_excel
#'
#' @export
ImportPTMs <- function(MaxQuant_table_name, samples_info = NULL, fasta_database = NA, prioritize_MaxQuant_names = TRUE, pattern_intensity = "Rep") {

  if (length(MaxQuant_table_name) != 1) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}
  if (is.na(MaxQuant_table_name)) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}
  if (!is.character(MaxQuant_table_name)) {stop("MaxQuant_table_name must be a character vector of length 1, indicating the name of the MaxQuant table files, in txt format")}
  if (!(endsWith(toupper(MaxQuant_table_name), ".TXT"))) {stop('MaxQuant_table_name must end with ".txt", and the file must be in the .txt format')}

  if (!is.null(samples_info)) {
    if (length(samples_info) != 1) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}

    if (!is.na(samples_info)) {
      if (!is.character(samples_info)) {stop("samples_info must be a character vector of length 1, indicating the name of the table in the current working directory containing information for each sample files")}
      if (!(endsWith(toupper(samples_info), ".TXT") | endsWith(toupper(samples_info), ".CSV") | endsWith(toupper(samples_info), ".XLSX"))) {stop('samples_info must end with ".txt", ".csv", or ".xlsx", and the file must be in that format')}
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



  proteinGroup_Raw <- read_tsv(MaxQuant_table_name, guess_max = Inf)


  if(!(all(c("Proteins", "Protein names", "Gene names", "id") %in% colnames(proteinGroup_Raw)))) {stop('MaxQuant_table must have at least the columns named: "Proteins", "Protein names", "Gene names", "id"')}

  if (!is.character(proteinGroup_Raw$`Proteins`)) stop('The column "Proteins" must contain character')

  if (!is.character(proteinGroup_Raw$`Protein names`)) stop('The column "Protein names" must contain character')

  if (!is.character(proteinGroup_Raw$`Gene names`)) stop('The column "Gene names" must contain character')


  if (!is.numeric(proteinGroup_Raw$id)) stop('The column "id" must contain numbers')

  column_intesities_names <- colnames(proteinGroup_Raw)[grepl(pattern_intensity, colnames(proteinGroup_Raw))]

  if (length(column_intesities_names)<1) {stop("no coulumns with the specified pattern_intensity!")}

  cat("\n___\n")
  cat(paste0("Considering the pattern_intensity '",pattern_intensity   ,"', the coulmns considered for the intensiteis are:\n '"))
  cat(paste0(column_intesities_names, collapse = "'\n '"))
  cat("'\n___\n")


  for (a in column_intesities_names) {
    if (!is.numeric(pull(proteinGroup_Raw, a))) {
      proteinGroup_Raw[,a] <- as.numeric(pull(proteinGroup_Raw, a))
    }
  }

  if ("protid" %in% colnames(proteinGroup_Raw)) {stop('Please, the MaxQuant table should not already have a column named "protid" as this function will create it, thanks!')}



  samplenames_raw_intensities <- column_intesities_names

  samplenames_raw_clean <- samplenames_raw_intensities

  if (any(duplicated(samplenames_raw_clean))) {stop("The names of samples must not be duplicated. Each of them must be unique!")}


  if (!is.null(samples_info)) {
    if (!is.na(samples_info)) {
      if (endsWith(toupper(samples_info), ".TXT")) {
        samples_info_loaded <- read_tsv(samples_info)
      } else if (endsWith(toupper(samples_info), ".CSV")) {
        samples_info_loaded <- read_csv(samples_info)
      } else if (endsWith(toupper(samples_info), ".XLSX")) {
        samples_info_loaded <- read_excel(samples_info, sheet = 1)
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
      samples_info_loaded <- tibble(Samples = samplenames_raw_clean)
    }
  } else {
    samples_info_loaded <- tibble(Samples = samplenames_raw_clean)
  }




  if (!is.null(fasta_database)) {
    if (!is.na(fasta_database)) {
      if (fasta_database == "human") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Human_ref.txt", package = "GetCoolProteopipe"))
      } else if (fasta_database == "mouse") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Mouse_ref.txt", package = "GetCoolProteopipe"))
      } else if (endsWith(fasta_database, ".txt")) {
        fasta_database_loaded <- read_tsv(fasta_database)
      } else if (endsWith(fasta_database, ".csv")) {
        fasta_database_loaded <- read_csv(fasta_database)
      }

      if (!(all(c("Accession", "Protein names", "Gene names") %in% colnames(fasta_database_loaded))))  {stop('The fasta database loaded must contain at least these columns: "Accession", "Protein names", and "Gene names"')}
      if (any(duplicated(pull(fasta_database_loaded, "Accession")))) {stop('The fasta database must not contain any duplicated in the "Accession" column')}
    }
  }


  if (any(duplicated(proteinGroup_Raw$id))) {

    id_deduplicated <- fix_duplicated(as.character(proteinGroup_Raw$id),
                                      zeros = TRUE,
                                      define_highest_for_zeros = NULL,
                                      start_with_zero = FALSE,
                                      exclude_the_first = FALSE,
                                      NA_as_character = TRUE)

  } else {
    id_deduplicated <- as.character(proteinGroup_Raw$id)
  }


  proteinGroup_Raw <- proteinGroup_Raw %>%
    add_column(protid = "prot",
               .before = 1) %>%
    mutate(protid = ifelse(id<10, paste0("prot00000", id_deduplicated),
                         ifelse(id<100, paste0("prot0000", id_deduplicated),
                                ifelse(id<1000, paste0("prot000", id_deduplicated),
                                       ifelse(id<10000, paste0(id_deduplicated),
                                              ifelse(id<100000, paste0("prot0", id_deduplicated),
                                                     paste0("prot", id_deduplicated)))))))




  #Add protein and gene names when missing if possible

  if ("Accession" %in% colnames(proteinGroup_Raw)) {warning('A column named "Accession" was already present in the imported table. Please not that now it has been completely replaced with each first code (before the ";") of "Proteins"')}

  proteinGroup_Raw_acc <-  mutate(proteinGroup_Raw, Accession = gsub(";.*", "", `Proteins`))

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


  #Keep only first name in protein names and gene names

  proteinGroup_Raw_added_cut <- mutate(proteinGroup_Raw_added,
                                       `Protein names` = gsub(";.*", "", `Protein names`),
                                       `Gene names` = gsub(";.*", "", `Gene names`))



  colnames_raw_intensities <- c("protid", samplenames_raw_intensities)

  colnames_LFQ_intensities <- colnames_raw_intensities

  colnames_others <- c("protid", colnames(proteinGroup_Raw_added_cut)[which(!(colnames(proteinGroup_Raw_added_cut)%in%colnames_raw_intensities[-1]) & !(colnames(proteinGroup_Raw_added_cut)=="protid"))])


  final_list <- list(quant_raw = select(proteinGroup_Raw_added_cut, all_of(colnames_raw_intensities)),
                     quant_LFQ = select(proteinGroup_Raw_added_cut, all_of(colnames_LFQ_intensities)),
                     proteinINFO = select(proteinGroup_Raw_added_cut, all_of(colnames_others)),
                     sampleINFO = samples_info_loaded)


  final_list$quant_LFQ[, colnames(final_list$quant_LFQ)[which(colnames(final_list$quant_LFQ) != "protid")]] <- as.numeric(NA)

  final_list$quant_raw <- mutate_at(final_list$quant_raw,
                                    colnames(final_list$quant_raw)[which(colnames(final_list$quant_raw) != "protid")],
                                    ~ replace(., is.nan(.), NA_real_))

  return(final_list)

}

