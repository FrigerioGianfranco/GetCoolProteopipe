#' Scale intensities
#'
#' In the intensities data table, it scale all the intensities as specified.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param subtract NULL, or "shift_median", or "mean", or "median". Specify it to subtract the mean or median value to all values. If "shift_median", an offset will be subtracted in order to center the medians to the global median.
#' @param divide NULL, or "sqrt_sd", or "sd", or "maxmin". Specify it to divide all values by the standard deviation, by the square root of the standard deviation, or by the range (max-min).
#' @param by_sample logical. If TRUE, the centering and scaling will be column-wise (per sample); if FALSE; row-wise (per protein).
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups. It needs to be passed only if subtract is "shift_median" and you if want to center each sample in a group to the global median of that group, instead of the overall global median.
#'
#' @return a GCPlist list with the scaled data intensity tables.
#'
#' @examples
#' \dontrun{
#'
#' GCPlist10 <- GCP_ScaleIntensities(GCPlist = GCPlist09,
#'                                   subtract = "shift_median",
#'                                   divide = "sqrt_sd",
#'                                   by_sample = TRUE)
#'
#'
#' }
#'
#'
#' @export
GCP_ScaleIntensities <- function(GCPlist, subtract = c("shift_median", "mean", "median"), divide = c("sqrt_sd", "sd", "maxmin"), by_sample = TRUE, name_column_groups = getOption("GetCoolProteopipe.name_column_groups")) {

  checkGCPlist(GCPlist)

  if (!is.null(subtract)) {

    if (!identical(tolower(subtract), c("shift_median", "mean", "median"))) {
      if (length(subtract) != 1) {stop('subtract must be one of "shift_median", "mean", "median"')}
      if (is.na(subtract)) {stop('subtract must be one of "shift_median", "mean", "median"')}
    }
    subtract <- tolower(subtract)
    subtract <- match.arg(subtract, c("shift_median", "mean", "median"))
  }

  if (!is.null(divide)) {

    if (!identical(tolower(divide), c("sqrt_sd", "sd", "maxmin"))) {
      if (length(divide) != 1) {stop('divide must be one of "sqrt_sd", "sd", "maxmin"')}
      if (is.na(divide)) {stop('divide must be one of "sqrt_sd", "sd", "maxmin"')}
    }
    divide <- tolower(divide)
    divide <- match.arg(divide, c("sqrt_sd", "sd", "maxmin"))
  }

  if (!is.logical(by_sample)) {stop("by_sample must be either TRUE or FALSE")}
  if (length(by_sample) != 1) {stop("by_sample must be either TRUE or FALSE")}
  if (is.na(by_sample)) {stop("by_sample must be either TRUE or FALSE")}

  if (subtract == "shift_median") {
    if (!is.null(name_column_groups)) {
      if (length(name_column_groups)!=1) {stop("name_column_groups must be NULL or a character of length 1")}
      if (!is.character(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1")}
      if (is.na(name_column_groups)) {stop("name_column_groups must be NULL or a character of length 1, not a NA")}
      cat(paste0("\n -- The name_column_groups considered is '", name_column_groups, "' --\n\n"))
      if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

      if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
        GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
      }
    } else {
      cat("\n -- The name_column_groups considered is NULL --\n\n")
    }
  }


  GCPoutput <- GCPlist

  if (by_sample) {

    if (!is.null(subtract)) {
      if (subtract == "shift_median") {
        if (is.null(name_column_groups)) {
          global_median <- median(map_dbl(GCPlist$intensities[, which(colnames(GCPlist$intensities)!="protid")], median))
        } else {
          global_medians_groups <- numeric(length = length(levels(pull(GCPlist$sampleINFO, name_column_groups))))
          names(global_medians_groups) <- levels(pull(GCPlist$sampleINFO, name_column_groups))
          for (grp in levels(pull(GCPlist$sampleINFO, name_column_groups))) {
            global_medians_groups[which(names(global_medians_groups)==grp)] <- median(map_dbl(GCPlist$intensities[, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)%in%pull(GCPlist$sampleINFO[which(pull(GCPlist$sampleINFO, name_column_groups)==grp),], 1))]], median))
          }
        }
      }
    }


    for (a in colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]) {

      if (is.null(subtract)) {
        subtract_by <- 0
      } else {
        if (subtract == "mean") {
          subtract_by <- mean(pull(GCPlist$intensities, a), na.rm = TRUE)
        } else if (subtract == "median") {
          subtract_by <- median(pull(GCPlist$intensities, a), na.rm = TRUE)
        } else if (subtract == "shift_median") {
          if (is.null(name_column_groups)) {
            subtract_by <- median(pull(GCPlist$intensities, a), na.rm = TRUE) - global_median
          } else {
            subtract_by <- median(pull(GCPlist$intensities, a), na.rm = TRUE) - global_medians_groups[as.character(pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1)==a)])]
          }
        }
      }


      if (is.null(divide)) {
        divide_by <- 1
      } else {
        if (divide == "sd") {
          divide_by <- sd(pull(GCPlist$intensities, a), na.rm = TRUE)
        } else if (divide == "sqrt_sd") {
          divide_by <- sqrt(sd(pull(GCPlist$intensities, a), na.rm = TRUE))
        } else if (divide == "maxmin") {
          divide_by <- max(pull(GCPlist$intensities, a), na.rm = TRUE) - min(pull(GCPlist$intensities, a), na.rm = TRUE)
        }
      }

      for (i in 1:length(pull(GCPlist$intensities, a))) {
        GCPoutput$intensities[i, a] <- (pull(GCPlist$intensities, a)[i] - subtract_by)/divide_by
      }
    }

  } else {

    for (i in 1:length(pull(GCPlist$intensities, 1))) {

      if (is.null(subtract)) {
        subtract_by <- 0
      } else {
        if (subtract == "mean") {
          subtract_by <- mean(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE)
        } else if (subtract == "median") {
          subtract_by <- median(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE)
        } else {
          stop("subtracting by shift_median will be computed only if by_sample is TRUE")
        }
      }


      if (is.null(divide)) {
        divide_by <- 1
      } else {
        if (divide == "sd") {
          divide_by <- sd(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE)
        } else if (divide == "sqrt_sd") {
          divide_by <- sqrt(sd(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE))
        } else if (divide == "maxmin") {
          divide_by <- max(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE) - min(as.vector(GCPlist$intensities[i, colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]], mode = "numeric"), na.rm = TRUE)
        }
      }


      for (a in colnames(GCPlist$intensities)[which(colnames(GCPlist$intensities)!= "protid")]) {

        GCPoutput$intensities[i, a] <- (pull(GCPlist$intensities, a)[i] - subtract_by)/divide_by

      }
    }
  }

  cat("\n______\nIntensities have been, ")
  if (by_sample) {
    cat("by sample")
  } else {
    cat("by protein")
  }
  if(!is.null(subtract)) {
    if (subtract == "mean") {
      cat(",\n subtracted by the mean")
    } else if (subtract == "median") {
      cat(",\n subtracted by the median")
    } else if (subtract == "shift_median") {
      if (is.null(name_column_groups)) {
        cat(",\n centered to the global median")
      } else {
        cat(paste0(",\n centered to the global median for each sample group (accorting to the '", name_column_groups, "' column)"))
      }
    }
  }
  if(!is.null(divide)) {
    if (divide == "sd") {
      cat(",\n divided by the standard deviation")
    } else if (divide == "sqrt_sd") {
      cat(",\n divided by the square root of the standard deviation")
    } else if (divide == "maxmin") {
      cat(",\n divided by the range (the difference betweein maximum and minimum values)")
    }
  }
  cat(".\n______\n")

  return(GCPoutput)
}



