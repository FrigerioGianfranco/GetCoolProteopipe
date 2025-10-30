#' Scale intensities
#'
#' In the quant_raw and the quantLFQ data intensity, it scale all the intensities as specified.
#'
#' @param GCPlist a list created with the ImportOutputMaxQuant function.
#' @param raw_or_LFQ one of the following: "raw", "LFQ". The scale will be performed only in the specified data intensities. Please note that it might not be a good idea to do so on LFQ data, so a warning is generated if LFQ is chosen.
#' @param subtract NULL, or "shift_median", or "mean", or "median". Specify it to subtract the mean or median value to all values. If "shift_median", an offset will be subtracted in order to center the medians to the global median.
#' @param divide NULL, or "sqrt_sd", or "sd", or "maxmin". Specify it to divide all values by the standard deviation, by the square root of the standard deviation, or by the range (max-min), for each protein.
#' @param by_sample logical. If TRUE, the centering and scaling will be column-wise (per sample); if FALSE; row-wise (per protein).
#' @param name_column_groups NULL or character of length 1. The name of the column of the sampleINFO table containing the sample groups. It needs to be passed only if subtract is "shift_median" and you want to center each sample in a group to the global median of that group, instead of the overall global median.
#'
#' @return a GCPlist list with the scaled data intensity tables.
#'
#' @export
GCP_ScaleIntensities <- function(GCPlist, raw_or_LFQ = getOption("GetCoolProteopipe.raw_or_LFQ"),
                                 subtract = c("shift_median", "mean", "median"), divide = c("sqrt_sd", "sd", "maxmin"), by_sample = TRUE, name_column_groups = NULL) {

  checkGCPlist(GCPlist)

  if (!identical(tolower(raw_or_LFQ), c("lfq", "raw"))) {
    if (length(raw_or_LFQ) != 1) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
    if (is.na(raw_or_LFQ)) {stop('raw_or_LFQ must be one of "raw", "LFQ"')}
  }
  raw_or_LFQ <- tolower(raw_or_LFQ)
  raw_or_LFQ <- match.arg(raw_or_LFQ, c("lfq", "raw"))

  if (raw_or_LFQ == "lfq") {
    cat("\n -- LFQ data are used --\n\n")
  } else if (raw_or_LFQ == "raw") {
    cat("\n -- raw data are used --\n\n")
  }

  if (raw_or_LFQ == "lfq") {
    warning("Are you sure is it a good idea to scale LFQ data?")
  }

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
      if (length(which(colnames(GCPlist$sampleINFO) == name_column_groups)) != 1) {stop("The name passed in name_column_groups must be a name of a column of the sampleINFO dataframe")}

      if (!is.factor(pull(GCPlist$sampleINFO, name_column_groups))) {
        GCPlist$sampleINFO[,name_column_groups] <- as.factor(pull(GCPlist$sampleINFO, name_column_groups))
      }
    }
  }



  GCPoutput <- GCPlist


  if (by_sample) {

    if (raw_or_LFQ == "raw") {
      if (!is.null(subtract)) {
        if (subtract == "shift_median") {
          if (is.null(name_column_groups)) {
            global_median <- median(map_dbl(GCPlist$quant_raw[, which(colnames(GCPlist$quant_raw)!="protid")], median))
          } else {
            global_medians_groups <- numeric(length = length(levels(pull(GCPlist$sampleINFO, name_column_groups))))
            names(global_medians_groups) <- levels(pull(GCPlist$sampleINFO, name_column_groups))
            for (grp in levels(pull(GCPlist$sampleINFO, name_column_groups))) {
              global_medians_groups[which(names(global_medians_groups)==grp)] <- median(map_dbl(GCPlist$quant_raw[, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)%in%pull(GCPlist$sampleINFO[which(pull(GCPlist$sampleINFO, name_column_groups)==grp),], 1))]], median))
            }
          }
        }
      }


      for (a in colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]) {

        if (is.null(subtract)) {
          subtract_by <- 0
        } else {
          if (subtract == "mean") {
            subtract_by <- mean(pull(GCPlist$quant_raw, a), na.rm = TRUE)
          } else if (subtract == "median") {
            subtract_by <- median(pull(GCPlist$quant_raw, a), na.rm = TRUE)
          } else if (subtract == "shift_median") {
            if (is.null(name_column_groups)) {
              subtract_by <- median(pull(GCPlist$quant_raw, a), na.rm = TRUE) - global_median
            } else {
              subtract_by <- median(pull(GCPlist$quant_raw, a), na.rm = TRUE) - global_medians_groups[as.character(pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1)==a)])]
            }
          }
        }


        if (is.null(divide)) {
          divide_by <- 1
        } else {
          if (divide == "sd") {
            divide_by <- sd(pull(GCPlist$quant_raw, a), na.rm = TRUE)
          } else if (divide == "sqrt_sd") {
            divide_by <- sqrt(sd(pull(GCPlist$quant_raw, a), na.rm = TRUE))
          } else if (divide == "maxmin") {
            divide_by <- max(pull(GCPlist$quant_raw, a), na.rm = TRUE) - min(pull(GCPlist$quant_raw, a), na.rm = TRUE)
          }
        }


        for (i in 1:length(pull(GCPlist$quant_raw, a))) {
          GCPoutput$quant_raw[i, a] <- (pull(GCPlist$quant_raw, a)[i] - subtract_by)/divide_by
        }
      }


    } else if (raw_or_LFQ == "lfq") {
      if (!is.null(subtract)) {
        if (subtract == "shift_median") {
          if (is.null(name_column_groups)) {
            global_median <- median(map_dbl(GCPlist$quant_LFQ[, which(colnames(GCPlist$quant_LFQ)!="protid")], median))
          } else {
            global_medians_groups <- numeric(length = length(levels(pull(GCPlist$sampleINFO, name_column_groups))))
            names(global_medians_groups) <- levels(pull(GCPlist$sampleINFO, name_column_groups))
            for (grp in levels(pull(GCPlist$sampleINFO, name_column_groups))) {
              global_medians_groups[which(names(global_medians_groups)==grp)] <- median(map_dbl(GCPlist$quant_LFQ[, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)%in%pull(GCPlist$sampleINFO[which(pull(GCPlist$sampleINFO, name_column_groups)==grp),], 1))]], median))
            }
          }
        }
      }

      for (a in colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]) {

        if (is.null(subtract)) {
          subtract_by <- 0
        } else {
          if (subtract == "mean") {
            subtract_by <- mean(pull(GCPlist$quant_LFQ, a), na.rm = TRUE)
          } else if (subtract == "median") {
            subtract_by <- median(pull(GCPlist$quant_LFQ, a), na.rm = TRUE)
          } else if (subtract == "shift_median") {
            if (is.null(name_column_groups)) {
              subtract_by <- median(pull(GCPlist$quant_LFQ, a), na.rm = TRUE) - global_median
            } else {
              subtract_by <- median(pull(GCPlist$quant_LFQ, a), na.rm = TRUE) - global_medians_groups[as.character(pull(GCPlist$sampleINFO, name_column_groups)[which(pull(GCPlist$sampleINFO, 1)==a)])]
            }
          }
        }


        if (is.null(divide)) {
          divide_by <- 1
        } else {
          if (divide == "sd") {
            divide_by <- sd(pull(GCPlist$quant_LFQ, a), na.rm = TRUE)
          } else if (divide == "sqrt_sd") {
            divide_by <- sqrt(sd(pull(GCPlist$quant_LFQ, a), na.rm = TRUE))
          } else if (divide == "maxmin") {
            divide_by <- max(pull(GCPlist$quant_LFQ, a), na.rm = TRUE) - min(pull(GCPlist$quant_LFQ, a), na.rm = TRUE)
          }
        }


        for (i in 1:length(pull(GCPlist$quant_LFQ, a))) {
          GCPoutput$quant_LFQ[i, a] <- (pull(GCPlist$quant_LFQ, a)[i] - subtract_by)/divide_by
        }
      }
    }





  } else {

    if (raw_or_LFQ == "raw") {

      for (i in 1:length(pull(GCPlist$quant_raw, 1))) {

        if (is.null(subtract)) {
          subtract_by <- 0
        } else {
          if (subtract == "mean") {
            subtract_by <- mean(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else if (subtract == "median") {
            subtract_by <- median(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else {
            stop("subtracting by shift_median will be computed only if by_sample is TRUE")
          }
        }


        if (is.null(divide)) {
          divide_by <- 1
        } else {
          if (divide == "sd") {
            divide_by <- sd(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else if (divide == "sqrt_sd") {
            divide_by <- sqrt(sd(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE))
          } else if (divide == "maxmin") {
            divide_by <- max(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE) - min(as.vector(GCPlist$quant_raw[i, colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          }
        }


        for (a in colnames(GCPlist$quant_raw)[which(colnames(GCPlist$quant_raw)!= "protid")]) {

          GCPoutput$quant_raw[i, a] <- (pull(GCPlist$quant_raw, a)[i] - subtract_by)/divide_by

        }
      }



    } else if (raw_or_LFQ == "lfq") {

      for (i in 1:length(pull(GCPlist$quant_LFQ, 1))) {

        if (is.null(subtract)) {
          subtract_by <- 0
        } else {
          if (subtract == "mean") {
            subtract_by <- mean(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else if (subtract == "median") {
            subtract_by <- median(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else {
            stop("subtracting by shift_median will be computed only if by_sample is TRUE")
          }
        }


        if (is.null(divide)) {
          divide_by <- 1
        } else {
          if (divide == "sd") {
            divide_by <- sd(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          } else if (divide == "sqrt_sd") {
            divide_by <- sqrt(sd(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE))
          } else if (divide == "maxmin") {
            divide_by <- max(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE) - min(as.vector(GCPlist$quant_LFQ[i, colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]], mode = "numeric"), na.rm = TRUE)
          }
        }


        for (a in colnames(GCPlist$quant_LFQ)[which(colnames(GCPlist$quant_LFQ)!= "protid")]) {

          GCPoutput$quant_LFQ[i, a] <- (pull(GCPlist$quant_LFQ, a)[i] - subtract_by)/divide_by

        }
      }
    }
  }


  return(GCPoutput)

}



