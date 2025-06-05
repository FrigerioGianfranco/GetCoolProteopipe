#' Importing the output of MaxQuant
#'
#' Starting from the output of the MaxQuant, it performs some cleaning such as filtering out proteins that are only identified by site, the reverse and the potential contaminants, then it adds protein names from a fasta database.
#'
#' @param MaxQuant_table_name a character vector of length 1 with the name of the MaxQuant table file in the current working directory, which must be in the .txt format.
#' @param samples_info NULL or NA or a character vector of length 1 with the name of the table in the current working directory, containing information for each sample. The table must be in txt, csv, or xslsx format. In particular, the first column of the table must contain the names of the samples exactly as they are in the MaxQuant table.
#' @param fasta_database NULL or NA or either "human" or "mouse". You can specify the fasta database to use to fill the missing Protein names and Gene names from the MaxQuant table
#' @param prioritize_MaxQuant_names logical. If TRUE and if a fasta_database is provided, the final "Protein names" and "Gene names" will be primarily taken from the  MaxQuant table (they will be taken from the fasta database only if missing). If FALSE, the opposite will happen.
#' @param remove_identified_by_site logical. Do you want to remove rows that contains "+" in the column "Only identified by site"?
#' @param remove_reverse logical. Do you want to remove rows that contains "+" in the column "Reverse"?
#' @param remove_potential_contaminant logical. Do you want to remove rows that contains "+" in the column "Potential contaminant"?
#'
#' @return a GCPlist, i.e.: a list with 4 dataframes (tibbles):#'
#'
#' - `quant_raw`: the raw intensities.
#'
#' - `quant_LFQ`: the LFQ intensities.
#'
#' - `proteinINFO`: all the information for each rows.
#'
#' - `sampleINFO`: the information about the samples.
#'
#'
#' @importFrom readxl read_excel
#'
#' @export
ImportOutputMaxQuant <- function(MaxQuant_table_name, samples_info = NULL, fasta_database = c(NA, "human", "mouse"), prioritize_MaxQuant_names = TRUE, remove_identified_by_site = TRUE, remove_reverse = TRUE, remove_potential_contaminant = TRUE) {

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
    if (!identical(tolower(fasta_database), c(NA, "human", "mouse"))) {
      if (length(fasta_database) != 1) {stop('fasta_database must be NULL or NA, or "human", or "mouse"')}
    }

    if (!is.na(fasta_database)) {
      if (tolower(fasta_database) == "na") {fasta_database <- NA}
    }

    fasta_database <- tolower(fasta_database)
    fasta_database <- match.arg(fasta_database, c(NA, "human", "mouse"))
  }


  if (!is.logical(prioritize_MaxQuant_names)) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}
  if (length(prioritize_MaxQuant_names) != 1) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}
  if (is.na(prioritize_MaxQuant_names)) {stop("prioritize_MaxQuant_names must be either TRUE or FALSE")}

  if (!is.logical(remove_identified_by_site)) {stop("remove_identified_by_site must be either TRUE or FALSE")}
  if (length(remove_identified_by_site) != 1) {stop("remove_identified_by_site must be either TRUE or FALSE")}
  if (is.na(remove_identified_by_site)) {stop("remove_identified_by_site must be either TRUE or FALSE")}

  if (!is.logical(remove_reverse)) {stop("remove_reverse must be either TRUE or FALSE")}
  if (length(remove_reverse) != 1) {stop("remove_reverse must be either TRUE or FALSE")}
  if (is.na(remove_reverse)) {stop("remove_reverse must be either TRUE or FALSE")}

  if (!is.logical(remove_potential_contaminant)) {stop("remove_potential_contaminant must be either TRUE or FALSE")}
  if (length(remove_potential_contaminant) != 1) {stop("remove_potential_contaminant must be either TRUE or FALSE")}
  if (is.na(remove_potential_contaminant)) {stop("remove_potential_contaminant must be either TRUE or FALSE")}




  ###1st object: proteinGroup_Raw
  #Load txt
  proteinGroup_Raw <- read_tsv(MaxQuant_table_name, guess_max = Inf)



  if(!(all(c("Protein IDs", "Protein names", "Gene names", "id") %in% colnames(proteinGroup_Raw)))) {stop('MaxQuant_table must have at least the columns named: "Protein IDs", "Protein names", "Gene names", "id"')}

  if (!is.character(proteinGroup_Raw$`Protein IDs`)) stop('The column "Protein IDs" must contain character')

  if (!is.character(proteinGroup_Raw$`Protein names`)) stop('The column "Protein names" must contain character')

  if (!is.character(proteinGroup_Raw$`Gene names`)) stop('The column "Gene names" must contain character')


  if (!is.numeric(proteinGroup_Raw$id)) stop('The column "id" must contain numbers')
  if (any(duplicated(proteinGroup_Raw$id))) {stop(paste0('The column "id" must not contain duplicates. In particular the following id are present more than once: ',
                                                         paste0(unique(proteinGroup_Raw$id[which(duplicated(proteinGroup_Raw$id))]), collapse = ", ")))}


  if (any(map_lgl(select(proteinGroup_Raw, all_of(colnames(proteinGroup_Raw)[which(grepl("intensity", tolower(colnames(proteinGroup_Raw))))])), ~any(is.na(.))))) {
    stop("all the Intensity columns must not contain missing values")
  }
  if (!(all(map_lgl(select(proteinGroup_Raw, all_of(colnames(proteinGroup_Raw)[which(grepl("intensity", tolower(colnames(proteinGroup_Raw))))])), is.numeric)))) {
    stop("all the Intensity columns must contain numerical data")
  }

  if ("protid" %in% colnames(proteinGroup_Raw)) {stop('Please, the MaxQuant table should not already have a column named "protid" as this function will create it, thanks!')}


  if ("intensity" %in% tolower(colnames(proteinGroup_Raw))) {
    colnames(proteinGroup_Raw)[which(tolower(colnames(proteinGroup_Raw))=="intensity")] <- "sum_raw_intensities"
  }

  samplenames_raw_intensities <- colnames(proteinGroup_Raw)[which(grepl("intensity", tolower(colnames(proteinGroup_Raw))) & !(grepl("lfq", tolower(colnames(proteinGroup_Raw)))))]

  if (length(samplenames_raw_intensities)==0) {warning("There are no raw intensities columns!")}

  samplenames_LFQ_intensities <- colnames(proteinGroup_Raw)[which(grepl("intensity", tolower(colnames(proteinGroup_Raw))) & grepl("lfq", tolower(colnames(proteinGroup_Raw))))]

  if (length(samplenames_LFQ_intensities)==0) {warning("There are no LFQ intensities columns!")}

  samplenames_raw_clean <- str_remove_all(samplenames_raw_intensities, regex("Intensity ", ignore_case = TRUE))

  samplenames_LFQ_clean <- str_remove_all(samplenames_LFQ_intensities, regex("LFQ Intensity ", ignore_case = TRUE))

  if (!(all(samplenames_LFQ_clean %in% samplenames_raw_clean) & all(samplenames_raw_clean %in% samplenames_LFQ_clean) & length(samplenames_raw_clean) == length(samplenames_LFQ_clean))) {
    stop('There are some differences among raw samples names and LFQ samples names...! \n Please, you might also also want to check if in the MaxQuant table the columns with raw intensities start with "Intensity ", and the columns with LFQ intensities start with "LFQ Intensity "')
  }

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


      if (all(pull(samples_info_loaded, 1) %in% samplenames_raw_intensities)) {
        samples_info_loaded[,1] <- str_remove_all(pull(samples_info_loaded, 1), regex("Intensity ", ignore_case = TRUE))
      } else if (all(pull(samples_info_loaded, 1) %in% samplenames_LFQ_intensities)) {
        samples_info_loaded[,1] <- str_remove_all(pull(samples_info_loaded, 1), regex("LFQ Intensity ", ignore_case = TRUE))
      } else if (!all(pull(samples_info_loaded, 1) %in% samplenames_raw_clean)) {
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
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Human_ref.txt", package = "GetCoolproteopipe"))

      } else if (fasta_database == "mouse") {
        fasta_database_loaded <- read_tsv(system.file("extdata", "Database_Mouse_ref.txt", package = "GetCoolproteopipe"))
      }

      if (!(all(c("Accession", "Protein names", "Gene names") %in% colnames(fasta_database_loaded))))  {stop('The fasta database loaded must contain at least these columns: "Accession", "Protein names", and "Gene names"')}
      if (any(duplicated(pull(fasta_database_loaded, "Accession")))) {stop('The fasta database must not contain any duplicated in the "Accession" column')}
    }
  }


  proteinGroup_Raw <- proteinGroup_Raw %>%
    add_column(protid = "prot",
               .before = 1) %>%
    mutate(protid = ifelse(id<10, paste0("prot00000", id),
                         ifelse(id<100, paste0("prot0000", id),
                                ifelse(id<1000, paste0("prot000", id),
                                       ifelse(id<10000, paste0("prot00", id),
                                              ifelse(id<100000, paste0("prot0", id),
                                                     paste0("prot", id)))))))


  cat("\n")
  cat("___\n")

  if (any(c(remove_identified_by_site, remove_reverse, remove_potential_contaminant))) {cat(paste0("At first, ", length(pull(proteinGroup_Raw, 1))  , " rows were imported\n"))}

  if (remove_identified_by_site) {
    if (!("Only identified by site" %in% colnames(proteinGroup_Raw))) {stop('"Only identified by site" must be a column of the MaxQuant_table')}
    if(any(proteinGroup_Raw$`Only identified by site`[which(!is.na(proteinGroup_Raw$`Only identified by site`))] != "+")) {stop('The column "Only identified by site" must contain only "+" or NA')}


    proteinGroup_Raw[which(is.na(pull(proteinGroup_Raw, "Only identified by site"))), "Only identified by site"] <- ""

    cat(paste0("...from which, ", sum(pull(proteinGroup_Raw, "Only identified by site") =="+"), " 'Only identified by site' were removed\n"))
    proteinGroup_Raw <-filter(proteinGroup_Raw, `Only identified by site` !="+")
  }

  if (remove_reverse) {
    if (!("Reverse" %in% colnames(proteinGroup_Raw))) {stop('"Reverse" must be a column of the MaxQuant_table')}
    if(any(proteinGroup_Raw$Reverse[which(!is.na(proteinGroup_Raw$Reverse))] != "+")) {stop('The column "Reverse" must contain only "+" or NA')}

    proteinGroup_Raw[which(is.na(pull(proteinGroup_Raw, "Reverse"))), "Reverse"] <- ""

    cat(paste0("...from which, ", sum(pull(proteinGroup_Raw, "Reverse") =="+"), " 'Reverse' were removed\n"))
    proteinGroup_Raw <-filter(proteinGroup_Raw, Reverse !="+")
  }

  if (remove_potential_contaminant) {
    if (!("Potential contaminant" %in% colnames(proteinGroup_Raw))) {stop('"Potential contaminant" must be a column of the MaxQuant_table')}
    if(any(proteinGroup_Raw$`Potential contaminant`[which(!is.na(proteinGroup_Raw$`Potential contaminant`))] != "+")) {stop('The column "Potential contaminant" must contain only "+" or NA')}

    proteinGroup_Raw[which(is.na(pull(proteinGroup_Raw, "Potential contaminant"))), "Potential contaminant"] <- ""

    cat(paste0("...from which, ", sum(pull(proteinGroup_Raw, "Potential contaminant") =="+"), " 'Potential contaminant' were removed\n"))
    proteinGroup_Raw <-filter(proteinGroup_Raw, `Potential contaminant` !="+")
  }

  if (any(c(remove_identified_by_site, remove_reverse, remove_potential_contaminant))) {cat(paste0("In the end, ", length(pull(proteinGroup_Raw, 1))  , " rows were kept\n"))} else {cat(paste0(length(pull(proteinGroup_Raw, 1))  , " rows were imported and kept\n"))}

  cat("___\n")

  #Add protein and gene names when missing if possible

  if ("Accession" %in% colnames(proteinGroup_Raw)) {warning('A column named "Accession" was already present in the imported table. Please not that now it has been completely replaced with each first code (before the ";") of "Protein IDs"')}

  proteinGroup_Raw_acc <-  mutate(proteinGroup_Raw, Accession = gsub(";.*", "", `Protein IDs`))

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

  colnames_LFQ_intensities <- c("protid", samplenames_LFQ_intensities)

  colnames_others <- c("protid", colnames(proteinGroup_Raw_added_cut)[which(!(colnames(proteinGroup_Raw_added_cut)%in%colnames_raw_intensities[-1]) & !(colnames(proteinGroup_Raw_added_cut)%in%colnames_LFQ_intensities[-1])  & !(colnames(proteinGroup_Raw_added_cut)=="protid"))])


  final_list <- list(quant_raw = select(proteinGroup_Raw_added_cut, all_of(colnames_raw_intensities)),
                     quant_LFQ = select(proteinGroup_Raw_added_cut, all_of(colnames_LFQ_intensities)),
                     proteinINFO = select(proteinGroup_Raw_added_cut, all_of(colnames_others)),
                     sampleINFO = samples_info_loaded)


  colnames(final_list$quant_raw)[-1] <- samplenames_raw_clean
  colnames(final_list$quant_LFQ)[-1] <- samplenames_LFQ_clean


  return(final_list)

}

