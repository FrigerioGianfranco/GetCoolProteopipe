#' Assign samples to groups.
#'
#' It assigns each samples to a desired groups.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param automatic_assignment One of the following: "no", "groupfirst", "replicatefirst". If "no" is selected, groups will be assigned considering the other following arguments; if "groupfirst" or "replicatefirst", the groups will be automatically created from sample names, considering the separator passed to separator_automatic_assignment, and an additional column called as name_column_replicates will be also created.
#' @param separator_automatic_assignment character of length 1. If automatic_assignment is "groupfirst" or "replicatefirst", this separator will be used to distinguish group from sample names. Example: if sample names are "control_01", "control_02", "disease_03", "disease_04"; automatic_assignment is "groupfirst"; and separator_automatic_assignment is "_", the samples "control_01", "control_02", will be assigned to the "control" group and the samples "disease_03", "disease_04" will be assigned to the "disease" group. An additional column, called as name_column_replicates, with "01", "02", "03", "04" will be also created.
#' @param sample_names character vector containing the existing name of samples (Will be considered if automatic_assignment is "no").
#' @param group_names character vector containing the groups corresponding to the samples passed to sample_names (Will be considered if automatic_assignment is "no").
#' @param sample_group_table alternatively, you can pass here a table with the sample names in the first column and the relative group in the second column. If you pass an argument here (and if automatic_assignment is "no"), this will be considered instead of sample_names and group_names.
#' @param name_column_groups character of length 1. The name of the new column that will contain the groups.
#' @param name_column_replicates character of length 1. The name of the new column that will contain the replicates (it will be used if automatic_assignment is not "no").
#' @param factorlevels NULL or character. You can specify here the levels of the groups (the reference group should be put as first element). If NULL, it will just apply the function as.factor.
#' @param controlgroup NULL or character. You can specify here the name of the control group(s) which will be put as the first level(s).
#'
#' @return the GCPlist with an additional column in the sampleINFO data frames containing the assigned groups.
#'
#'
#' @examples
#' \dontrun{
#'
#' GCPlist02 <- GCP_AssignGroups(GCPlist = GCPlist01,
#'                               automatic_assignment = "groupfirst",
#'                               separator_automatic_assignment = "_",
#'                               name_column_groups = "Condition",
#'                               controlgroup = "S")
#'
#'
#' ## to just modify the order of levels of a group already present:
#'
#' GCPlist02bis <- GCP_AssignGroups(GCPlist = GCPlist02,
#'                                  group_names = GCPlist02$sampleINFO$Condition,
#'                                  name_column_groups = "Condition",
#'                                  factorlevels = c("S", "V"))
#'
#' }
#'
#'
#' @export
GCP_AssignGroups <- function(GCPlist, automatic_assignment = c("no", "groupfirst", "replicatefirst"), separator_automatic_assignment = "_", sample_names = pull(GCPlist$sampleINFO, 1), group_names = rep("not_assigned", length(sample_names)), sample_group_table = NULL, name_column_groups = "Condition", name_column_replicates = "Replicate", factorlevels = NULL, controlgroup = NULL) {

  checkGCPlist(GCPlist)

  if (!identical(tolower(automatic_assignment), c("no", "groupfirst", "replicatefirst"))) {
    if (length(automatic_assignment) != 1) {stop('automatic_assignment must be one of "no", "groupfirst", "replicatefirst"')}
    if (is.na(automatic_assignment)) {stop('automatic_assignment must be one of "no", "groupfirst", "replicatefirst"')}
  }
  automatic_assignment <- tolower(automatic_assignment)
  automatic_assignment <- match.arg(automatic_assignment, c("no", "groupfirst", "replicatefirst"))

  if (length(name_column_groups)!=1) {stop("name_column_groups must be a character of length 1")}
  if (!is.character(name_column_groups)) {stop("name_column_groups must be a character of length 1")}
  if (is.na(name_column_groups)) {stop("name_column_groups must a character of length 1, not a missing value")}

  if (!is.null(factorlevels)) {
    if (!is.character(factorlevels)) {stop("if not NULL, factorlevels must be a character")}
    if (length(factorlevels) < 1) {stop("if not NULL, factorlevels must not be empty!")}
    if (any(is.na(factorlevels))) {stop("if not NULL, factorlevels must be a character with no missing values")}
  }

  if (!is.null(controlgroup)) {
    if (!is.character(controlgroup)) {stop("if not NULL, controlgroup must be a character")}
    if (length(controlgroup) < 1) {stop("if not NULL, controlgroup must not be empty!")}
    if (any(is.na(controlgroup))) {stop("if not NULL, controlgroup must be a character with no missing values")}
  }


  if (automatic_assignment != "no") {
    if (length(separator_automatic_assignment) != 1) {stop("separator_automatic_assignment must be a character of length 1 with the separator to separe group names from sample names")}
    if (!is.character(separator_automatic_assignment)) {stop("separator_automatic_assignment must be a character of length 1 with the separator to separe group names from sample names")}
    if (is.na(separator_automatic_assignment)) {stop("separator_automatic_assignment must be a character of length 1 with the separator to separe group names from sample names, not a NA")}

    if (length(name_column_replicates)!=1) {stop("name_column_replicates must be a character of length 1")}
    if (!is.character(name_column_replicates)) {stop("name_column_replicates must be a character of length 1")}
    if (is.na(name_column_replicates)) {stop("name_column_replicates must a character of length 1, not a missing value")}

    names_first_parts <- str_remove(pull(GCPlist[["sampleINFO"]], 1), paste0(separator_automatic_assignment, ".*"))
    names_second_parts <- str_remove(pull(GCPlist[["sampleINFO"]], 1), paste0(names_first_parts, separator_automatic_assignment))

    GCPoutput <- GCPlist

    GCPoutput$sampleINFO[, name_column_groups] <- rep(as.character(NA), length(pull(GCPoutput$sampleINFO, 1)))
    GCPoutput$sampleINFO[, name_column_replicates] <- rep(as.character(NA), length(pull(GCPoutput$sampleINFO, 1)))

    if (automatic_assignment == "groupfirst") {
      GCPoutput$sampleINFO[, name_column_groups] <- names_first_parts
      GCPoutput$sampleINFO[, name_column_replicates] <- names_second_parts
    }
    if (automatic_assignment == "replicatefirst") {
      GCPoutput$sampleINFO[, name_column_groups]  <- names_second_parts
      GCPoutput$sampleINFO[, name_column_replicates] <- names_first_parts
    }

    if (any(duplicated(pull(GCPoutput$sampleINFO, 1)))) {
      stop(paste0('Duplicated names have been introduced in samples names. In particular, the following are duplicated: "', paste0(unique(pull(GCPoutput$sampleINFO, 1)[which(duplicated(pull(GCPoutput$sampleINFO, 1)))]), collapse = '", "'), '"'))
    }
  }

  if (automatic_assignment == "no") {
    if (is.null(sample_group_table)) {

      if (any(is.na(sample_names))) {stop("sample_names must not contain missing values")}

      if (length(sample_names) != length(group_names)) {stop("old_names and new_names must have the same length")}


      sample_names <- as.character(sample_names)
      group_names <- as.character(group_names)

      if (!all(sample_names %in% pull(GCPlist$sampleINFO, 1))) {
        stop(paste0("Names passed to sample_names must be names of samples. The following sample names you passed are not present:\n",
                    paste0(sample_names[which(!sample_names %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
      }

      if (any(duplicated(sample_names))) {
        stop(paste0("sample_names must not contain duplicated. The following sample names you passed are duplicated:\n",
                    paste0(unique(sample_names[which(duplicated(sample_names))]), collapse = "\n")))
      }

      if (any(is.na(group_names))) {
        group_names[which(is.na(group_names))] <- "not_assigned"
      }


      the_sample_group_table <- tibble(Sample_names = sample_names,
                                       Group_names = group_names)

      if (is.null(name_column_groups)) {name_column_groups <- deparse(substitute(group_names))}

    } else {

      sample_group_table <- as.data.frame(sample_group_table)

      if (ncol(sample_group_table) != 2) {stop("sample_group_table must contain exactly two columns")}

      if (any(is.na(pull(sample_group_table, 1)))) {stop("the first column of sample_group_table must not contain missing values")}

      the_sample_group_table <- tibble(Sample_names = as.character(pull(sample_group_table, 1)),
                                       Group_names = as.character(pull(sample_group_table, 2)))

      if (!all(the_sample_group_table$Sample_names %in% pull(GCPlist$sampleINFO, 1))) {
        stop(paste0("Names passed to the first column of sample_group_table must be names of samples. The following sample names you passed are not present:\n",
                    paste0(the_sample_group_table$Sample_names[which(!the_sample_group_table$Sample_names %in% pull(GCPlist$sampleINFO, 1))], collapse = "\n")))
      }

      if (any(duplicated(the_sample_group_table$Sample_names))) {
        stop(paste0("Names passed to the first column of sample_group_table must not contain duplicated. The following sample names you passed are duplicated:\n",
                    paste0(unique(the_sample_group_table$Sample_names[which(duplicated(the_sample_group_table$Sample_names))]), collapse = "\n")))
      }

      if (any(is.na(the_sample_group_table$Group_names))) {
        the_sample_group_table$Group_names[which(is.na(the_sample_group_table$Group_names))] <- "not_assigned"


        if (is.null(name_column_groups)) {name_column_groups <- colnames(sample_group_table)[2]}
      }
    }


    if (name_column_groups %in% colnames(GCPlist$sampleINFO)) {
      warning(paste0("The sampleINFO data frame already contaied a column named ", name_column_groups, ". Please, note that it has now been completely replaced. Consider passing another name to name_column_groups if you don't want to lose it."))
    }


    GCPoutput <- GCPlist

    GCPoutput$sampleINFO[, name_column_groups] <- rep("not_assigned", nrow(GCPoutput$sampleINFO))

    for (i in 1:nrow(GCPoutput$sampleINFO)) {

      if (pull(GCPoutput$sampleINFO, 1)[i] %in% the_sample_group_table$Sample_names) {
        GCPoutput$sampleINFO[i, name_column_groups] <- the_sample_group_table$Group_names[which(the_sample_group_table$Sample_names == pull(GCPoutput$sampleINFO, 1)[i])]
      }
    }
  }




  if (is.null(factorlevels)) {
    GCPoutput$sampleINFO[, name_column_groups] <- as.factor(pull(GCPoutput$sampleINFO, name_column_groups))
  } else {
    GCPoutput$sampleINFO[, name_column_groups] <- factor(pull(GCPoutput$sampleINFO, name_column_groups), levels = factorlevels)
    if (any(is.na(pull(GCPoutput$sampleINFO, name_column_groups)))) {
      warning("Missing values introduced in the assigned groups. Maybe you need to check the factorlevels argument!")
    }
  }

  if (!is.null(controlgroup)) {

    if (!all(controlgroup%in%levels(pull(GCPoutput$sampleINFO, name_column_groups)))) {warning(paste0('The following names you passed in controlgroup were not considered as they are not valid names of the groups:\n ', paste0(controlgroup[which(!controlgroup%in%levels(pull(GCPoutput$sampleINFO, name_column_groups)))], collapse = " ")))}

    controlgroup <- controlgroup[which(controlgroup%in%levels(pull(GCPoutput$sampleINFO, name_column_groups)))]

    if (length(controlgroup) > 0) {

      GCPoutput$sampleINFO[, name_column_groups] <- factor(as.character(pull(GCPoutput$sampleINFO, name_column_groups)),
                                                           levels = c(controlgroup, levels(pull(GCPoutput$sampleINFO, name_column_groups))[which(levels(pull(GCPoutput$sampleINFO, name_column_groups)) != controlgroup)]))
    }
  }



  cat("\n")
  cat("______\n")
  cat(paste0("The groups have been assigned in the sampleINFO dataframes in the column named ", name_column_groups, ", in this way:\n"))

  if (automatic_assignment == "no") {
    print(as.data.frame(GCPoutput$sampleINFO[, c(1, which(colnames(GCPoutput$sampleINFO) == name_column_groups))]))
  } else {
    print(as.data.frame(GCPoutput$sampleINFO[, c(1, which(colnames(GCPoutput$sampleINFO) == name_column_groups), which(colnames(GCPoutput$sampleINFO) == name_column_replicates))]))
  }


  cat(paste0("\n\nThe order of the levels is the following: ", paste0(levels(pull(GCPoutput$sampleINFO, name_column_groups)), collapse = " ")))


  cat("\n______\n")

  return(GCPoutput)

}
