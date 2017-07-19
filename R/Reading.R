globalVariables(c("tmpdata"))
ReadInputs <- function (varnames, climatedata, constants, stopmissing, 
                        timestep = "daily", 
                        interp_missing_days = FALSE, 
                        interp_missing_entries = FALSE, 
                        interp_abnormal = FALSE,
                        missing_method = NULL, abnormal_method = NULL) 
{
  # compile date variable 
  if ("Year" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Year'")
  }
  if ("Month" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Month'")
  }
  if ("Day" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Day'")
  }
  if (timestep == "subdaily") {
    Date.subdaily <- strptime(paste(climatedata$Day, "/", 
                                    climatedata$Month, "/", climatedata$Year, " ", climatedata$Hour, 
                                    sep = ""), "%d/%m/%Y %H")
    Date.daily <- unique(as.Date(Date.subdaily, "%d/%m/%y"))
    Date.daily <- Date.daily[which(!is.na(Date.daily))]
    Date.monthly <- unique(as.yearmon(Date.subdaily, "%d/%m/%y"))
    J.temp <- zoo(as.numeric(format(Date.daily,"%j")), as.Date(Date.daily))
    i.temp <- unique(Date.monthly)
    i <- (i.temp - trunc(i.temp)) * 12 + 1
    Ndays.temp <- zoo(as.numeric(format(Date.daily,"%d")), as.Date(Date.daily))
    dateagg <- Date.subdaily
    StdDate.daily <- seq.Date(as.Date(Date.daily[1]), as.Date(Date.daily[length(Date.daily)]), 
                              by = "day")
    Missing_DateIndex.daily <- as.Date(setdiff(StdDate.daily, 
                                               as.Date(Date.daily)))
  } else if (timestep == "daily") {
    Date.daily <- strptime(paste(climatedata$Day, "/", climatedata$Month, 
                                 "/", climatedata$Year, sep = ""), "%d/%m/%Y")
    Date.monthly <- unique(as.yearmon(Date.daily, "%d/%m/%y"))
    J.temp <- zoo(as.numeric(format(Date.daily,"%j")), as.Date(Date.daily))
    i.temp <- unique(Date.monthly)
    i <- (i.temp - trunc(i.temp)) * 12 + 1
    Ndays.temp <- zoo(as.numeric(format(Date.daily,"%d")), as.Date(Date.daily))
    dateagg <- as.Date(Date.daily)
    StdDate.daily <- seq.Date(as.Date(Date.daily[1]), as.Date(Date.daily[length(Date.daily)]), 
                              by = "day")
    Missing_DateIndex.daily <- as.Date(setdiff(StdDate.daily, 
                                               as.Date(Date.daily)))
  }
  Stdzoo <- zoo(StdDate.daily, as.Date(StdDate.daily))
  if (length(Missing_DateIndex.daily) > 0) {
    message(paste("Warning: Number of missing date indices: ", 
                  length(Missing_DateIndex.daily), " days", sep = ""))
    message(paste("% missing date indices: ", signif(length(Missing_DateIndex.daily)/length(StdDate.daily), 
                                                     digits = -3), "%", sep = ""))
    if (length(Missing_DateIndex.daily) >= stopmissing[1]/100 * 
        nrow(climatedata)) {
      stop("missing date indices exceeds ", stopmissing[1], 
           "%, please use high quality data for calculation")
    }
    
    if (interp_missing_days == T) {
      
      Date.daily <- StdDate.daily
      J <- merge(J.temp, Stdzoo, all = TRUE, 
                      fill = NA)$J.temp
      J[which(is.na(J))] <- as.numeric(format(as.Date(time(J[which(is.na(J))])),"%j"))
      Ndays <- merge(Ndays.temp, Stdzoo, all = TRUE, 
                          fill = NA)$Ndays.temp
      Ndays.temp[which(is.na(Ndays.temp))] <- as.numeric(format(as.Date(time(Ndays.temp[which(is.na(Ndays.temp))])),"%d"))

      message(paste("All climate variables for missing dates will be interpolated with ", 
                    missing_method, sep = ""))
    } else {
      J = J.temp
      message("NA will be filled in for all climate variables for missing dates")
    }
  } else {
     J = J.temp
  }
  
  #dateagg <- as.Date(StdDate.daily)
  Ndays <- aggregate(Ndays.temp, as.yearmon, FUN = max)
  if (is.na(as.numeric(stopmissing[1])) | is.na(as.numeric(stopmissing[2])) | 
      is.na(as.numeric(stopmissing[3]))) {
    message("Please use three numeric values for the maximum allowable percentages of: ")
    message("1. missing date indices to the total number of days")
    message("2. missing data entries to the total number of data entries for each climate variable")
    stop("3. continuous missing data entries to the total number of data entries for each climate variable")
  } else {
    if (length(stopmissing) != 3) {
      stop("Please input a vector of length 3 for argument 'stopmissing'")
    } else {
      for (counter in 1:2) {
        if (as.numeric(stopmissing[counter]) < 1 | as.numeric(stopmissing[counter]) > 
            99) {
          stop("Please use values between 1 and 99 for the maximum allowable percentage of date indices/missing data entries")
        }
      }
    }
  }
  message(paste("The maximum acceptable percentage of date indices is", 
                stopmissing[1], "%"))
  message(paste("The maximum acceptable percentage of missing data is", 
                stopmissing[2], "%"))
  message(paste("The maximum acceptable percentage of continuous missing data is", 
                stopmissing[3], "%"))
  
  Tmax = NULL
  Tmin = NULL
  RHmax = NULL
  RHmin = NULL
  Rs = NULL
  u2 = NULL
  uz = NULL
  
  # Read in Tmax
  if ("Tmax" %in% varnames) {
    if (timestep != "daily") {
      stop("Variable Tmax can only be read in when date is recorded at daily timescale")
    }
    
    if ("Tmax" %in% (colnames(climatedata))) {
      Tmax.temp <- zoo(as.vector(climatedata$Tmax), 
                       dateagg)
      if ("TRUE" %in% (is.na(climatedata$Tmax))) {
        message("Warning: missing values in 'Tmax' (daily maximum temperature)")
        message(paste("Number of missing values in Tmax: ", 
                      sum(is.na(climatedata$Tmax))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmax))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Tmax))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Tmax)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
          stop("missing data of Tmax exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Tmax exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Tmax.temp <- ReadInput_InterpMissing("Tmax.temp", 
                                               Tmax.temp, missing_method)
          if (is.null(Tmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(Tmax.temp) > 100|as.vector(Tmax.temp) < (-50)  )) > 
          0) {
        message(paste("Number of day increments when Tmax has errors (Tmax > 100 deg or < -50 deg): ", 
                      length(which(as.vector(Tmax.temp) > 100 |as.vector(Tmax.temp) < (-50)))))
        if (interp_abnormal == T) {
          Tmax.temp <- ReadInput_InterpAbnormal("Tmax.temp", 
                                                upperdata = NULL, Tmax.temp, abnormal_method)
          if (is.null(Tmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(Missing_DateIndex.daily) > 0) {
        Tmax.temp <- merge(Tmax.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmax.temp
        if (interp_missing_days == T) {
          Tmax <- ReadInput_InterpMissing("Tmax.temp", 
                                          Tmax.temp, missing_method)
          if (is.null(Tmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmax <- Tmax.temp
        }
      } else {
        Tmax <- Tmax.temp
      }
    } else {
      stop("Missing variable of Tmax in input data")
    }
  }
  # Read in Tmin
  if ("Tmin" %in% varnames) {
    if (timestep != "daily") {
      stop("Variable Tmin can only be read in when date is recorded at daily timescale")
    }
    
    if ("Tmin" %in% (colnames(climatedata))) {
      Tmin.temp <- zoo(as.vector(climatedata$Tmin), 
                       dateagg)
      if ("TRUE" %in% (is.na(climatedata$Tmin))) {
        message("Warning: missing values in 'Tmin' (daily minimum temperature)")
        message(paste("Number of missing values in Tmin: ", 
                      sum(is.na(climatedata$Tmin))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmin))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Tmin))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$Tmin)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
          stop("missing data of Tmin exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in Tmin exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          Tmin.temp <- ReadInput_InterpMissing("Tmin.temp", 
                                               Tmin.temp, missing_method)
          if (is.null(Tmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      
      if (length(which(as.vector(Tmin.temp-Tmax)>0 | as.vector(Tmin.temp)<(-50) )) > 0) {
        message(paste("Number of day increments when Tmin has errors (Tmin > Tmax or Tmin < -50 deg): ", 
                      length(which(as.vector(Tmin.temp-Tmax)>0 | as.vector(Tmin.temp)<(-50)  ))))
        if (interp_abnormal == T) {
          Tmin.temp <- ReadInput_InterpAbnormal("Tmin.temp", 
                                                upperdata = Tmax, Tmin.temp, 
                                                abnormal_method)
          if (is.null(Tmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      
      if (length(Missing_DateIndex.daily) > 0) {
        Tmin.temp <- merge(Tmin.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Tmin.temp
        if (interp_missing_days == T) {
          Tmin <- ReadInput_InterpMissing("Tmin.temp", 
                                          Tmin.temp, missing_method)
          if (is.null(Tmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          Tmin <- Tmin.temp
        }
      } else {
        Tmin <- Tmin.temp
      }
    } else {
      stop("Missing variable of Tmin in input data")
    }
  }
  # Read in Temp (subdaily)
  if (!"Tmax" %in% varnames & !"Tmin" %in% varnames ) {
    if ("Temp" %in% varnames) {
      if (timestep != "subdaily") {
        stop("Variable Temp can only be read in when date is recorded at subdaily timescale")
      }
      if ("Temp" %in% (colnames(climatedata))) {
        temp.temp <- zoo(as.vector(climatedata$Temp), 
                         dateagg)
        
        if ("TRUE" %in% (is.na(climatedata$Temp))) {
          message("Warning: missing values in 'Temp' (sub-daily temperature)")
          message(paste("Number of missing values in Temp: ", 
                        sum(is.na(climatedata$Temp))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Temp))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Temp))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Temp)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
            stop("missing data of Temp exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Temp exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            temp.temp <- ReadInput_InterpMissing("temp.temp", 
                                                 temp.temp, missing_method)
            if (is.null(temp.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(temp.temp) > 100|as.vector(temp.temp) < (-50)  )) > 
            0) {
          message(paste("Number of data entries where Temp has errors (Temp > 100 deg or < -50 deg): ", 
                        length(which(as.vector(temp.temp) > 100|as.vector(temp.temp) < (-50)))))
          if (interp_abnormal == T) {
            temp.temp <- ReadInput_InterpAbnormal("temp.temp", 
                                                  upperdata = NULL, temp.temp, abnormal_method)
            if (is.null(temp.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        Tmin.temp <- aggregate(temp.temp, as.Date(dateagg, 
                                                  "%d/%m/%y"), FUN = min)
        if (length(Missing_DateIndex.daily) > 0) {
          Tmin.temp <- merge(Tmin.temp, Stdzoo, all = TRUE, 
                             fill = NA)$Tmin.temp
          if (interp_missing_days == T) {
            Tmin <- ReadInput_InterpMissing("Tmin.temp", 
                                            Tmin.temp,  missing_method)
            if (is.null(Tmin)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Tmin <- Tmin.temp
          }
        } else {
          Tmin <- Tmin.temp
        }
        
        Tmax.temp <- aggregate(temp.temp, as.Date(dateagg, 
                                                  "%d/%m/%y"), FUN = max)
        if (length(Missing_DateIndex.daily) > 0) {
          Tmax.temp <- merge(Tmax.temp, Stdzoo, all = TRUE, 
                             fill = NA)$Tmax.temp
          if (interp_missing_days == T) {
            Tmax <- ReadInput_InterpMissing("Tmax.temp", 
                                            Tmax.temp,  missing_method)
            if (is.null(Tmax)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Tmax <- Tmax.temp
          }
        } else {
          Tmax <- Tmax.temp
        }
      } else {
        stop("Missing variable of Temp in input data")
      }
      
    } else {
    stop("Must have at least one of either daily 'Tmax' and 'Tmin', or subdaily 'Temp' for estimating ET")
    }
  }
  
  # Read in Tdew
  if ("Tdew" %in% varnames) {
    if (timestep == "daily") {
      if ("Tdew" %in% (colnames(climatedata))) {
        Tdew.temp <- zoo(as.vector(climatedata$Tdew), 
                         dateagg)
        if ("TRUE" %in% (is.na(climatedata$Tdew))) {
          message("Warning: missing values in 'Tdew' (daily dew point temperature)")
          message(paste("Number of missing values in Tdew: ", 
                        sum(is.na(climatedata$Tdew))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Tdew))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Tdew))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Tdew)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Tdew exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Tdew exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Tdew.temp <- ReadInput_InterpMissing("Tdew.temp", 
                                                 Tdew.temp, missing_method)
            if (is.null(Tdew.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        
        if (length(which(as.vector(Tdew.temp) > 100 |as.vector(Tdew.temp)<(-50) )) > 0) {
          message(paste("Number of day increments when Tdew has errors (Tdew > 100 deg and Tdew < -50 deg): ", 
                        length(which(as.vector(Tdew.temp) > 100 | as.vector(Tdew.temp)<(-50)  ))))
          if (interp_abnormal == T) {
            Tdew.temp <- ReadInput_InterpAbnormal("Tdew.temp", 
                                                  upperdata = NULL, Tdew.temp, 
                                                  abnormal_method)
            if (is.null(Tdew.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        
        if (length(Missing_DateIndex.daily) > 0) {
          Tdew.temp <- merge(Tdew.temp, Stdzoo, all = TRUE, 
                             fill = NA)$Tdew.temp
          if (interp_missing_days == T) {
            Tdew <- ReadInput_InterpMissing("Tdew.temp", 
                                            Tdew.temp, missing_method)
            if (is.null(Tdew)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Tdew <- Tdew.temp
          }
        } else {
          Tdew <- Tdew.temp
        }
      } else {
        stop("Missing variable of Tdew in input data")
      }
    } else if (timestep == "subdaily") {
      if ("Tdew" %in% (colnames(climatedata))) {
        Tdew.temp <- zoo(as.vector(climatedata$Tdew), 
                         dateagg)
        if ("TRUE" %in% (is.na(climatedata$Tdew))) {
          message("Warning: missing values in 'Tdew' (sub-daily dew point temperature)")
          message(paste("Number of missing values in Tdew: ", 
                        sum(is.na(climatedata$Tdew))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Tdew))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Tdew))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Tdew)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Tdew exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Tdew exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Tdew.temp <- ReadInput_InterpMissing("Tdew.temp", 
                                                 Tdew.temp, missing_method)
            if (is.null(Tdew.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        
        if (length(which(as.vector(Tdew.temp) > 100 |as.vector(Tdew.temp)<(-50) )) > 0) {
          message(paste("Number of day increments when Tdew has errors (Tdew > 100 deg and Tdew < -50 deg): ", 
                        length(which(as.vector(Tdew.temp) > 100 | as.vector(Tdew.temp)<(-50)  ))))
          if (interp_abnormal == T) {
            Tdew.temp <- ReadInput_InterpAbnormal("Tdew.temp", 
                                                  upperdata = NULL, Tdew.temp, 
                                                  abnormal_method)
            if (is.null(Tdew.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        Tdew.temp <- aggregate(Tdew.temp, as.Date(dateagg, 
                                                  "%d/%m/%y"), FUN = mean)
        if (length(Missing_DateIndex.daily) > 0) {
          Tdew.temp <- merge(Tdew.temp, Stdzoo, all = TRUE, 
                             fill = NA)$Tdew.temp
          if (interp_missing_days == T) {
            Tdew <- ReadInput_InterpMissing("Tdew.temp", 
                                            Tdew.temp, missing_method)
            if (is.null(Tdew)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Tdew <- Tdew.temp
          }
        } else {
          Tdew <- Tdew.temp
        }
      } else {
        stop("Missing variable of Tdew in input data")
      }
    }
    
    
  }
  
  # Read in RHmax
  if ("RHmax" %in% varnames) {
    if (timestep != "daily") {
      stop("Variable RHmax can only be read in when date is recorded at daily timescale")
    }
    
    if ("RHmax" %in% (colnames(climatedata))) {
      RHmax.temp <- zoo(as.vector(climatedata$RHmax), 
                        dateagg)
      if ("TRUE" %in% (is.na(climatedata$RHmax))) {
        message("Warning: missing values in 'RHmax' (daily maximum relative humidity)")
        message(paste("Number of missing values in RHmax: ", 
                      sum(is.na(climatedata$RHmax))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmax))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RHmax))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RHmax)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
          stop("missing data of RHmax exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RHmax exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RHmax.temp <- ReadInput_InterpMissing("RHmax.temp", 
                                                RHmax.temp, missing_method)
          if (is.null(RHmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(which(as.vector(RHmax.temp) > 100|as.vector(RHmax.temp) < 0  )) > 
          0) {
        message(paste("Number of day increments when RHmax has errors (RHmax > 100% deg or < 0%): ", 
                      length(which(as.vector(RHmax.temp) > 100 |as.vector(RHmax.temp) < 0))))
        if (interp_abnormal == T) {
          RHmax.temp <- ReadInput_InterpAbnormal("RHmax.temp", 
                                                 upperdata = NULL, RHmax.temp, abnormal_method)
          if (is.null(RHmax.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      if (length(Missing_DateIndex.daily) > 0) {
        RHmax.temp <- merge(RHmax.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmax.temp
        if (interp_missing_days == T) {
          RHmax <- ReadInput_InterpMissing("RHmax.temp", 
                                           RHmax.temp, missing_method)
          if (is.null(RHmax)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmax <- RHmax.temp
        }
      } else {
        RHmax <- RHmax.temp
      }
    } else {
      stop("Missing variable of RHmax in input data")
    }
  }
  # Read in RHmin
  if ("RHmin" %in% varnames) {
    if (timestep != "daily") {
      stop("Variable RHmin can only be read in when date is recorded at daily timescale")
    }
    
    if ("RHmin" %in% (colnames(climatedata))) {
      RHmin.temp <- zoo(as.vector(climatedata$RHmin), 
                        dateagg)
      if ("TRUE" %in% (is.na(climatedata$RHmin))) {
        message("Warning: missing values in 'RHmin' (daily mainimum relative humidity)")
        message(paste("Number of missing values in RHmin: ", 
                      sum(is.na(climatedata$RHmin))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmin))/nrow(climatedata) * 
                                                   100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RHmin))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for (counter in 2:nrow(df)) {
          df$zcount[counter] <- ifelse(df$x[counter] == 
                                         0, df$zcount[counter - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", 
                      signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                      "%"))
        if (sum(is.na(climatedata$RHmin)) >= stopmissing[2]/100 * 
            nrow(climatedata)) {
          stop("missing data of RHmin exceeds ", 
               stopmissing[2], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
            stop("Maximum duration of missing data in RHmin exceeds ", 
                 stopmissing[3], "% of total data duration, please use high quality data for calculation")
          }
        }
        if (interp_missing_entries == T) {
          RHmin.temp <- ReadInput_InterpMissing("RHmin.temp", 
                                                RHmin.temp, missing_method)
          if (is.null(RHmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      
      if (length(which(as.vector(RHmin.temp-RHmax)>0 | as.vector(RHmin.temp)<0 )) > 0) {
        message(paste("Number of day increments when RHmin has errors (RHmin > RHmax or RHmin < 0%): ", 
                      length(which(as.vector(RHmin.temp-RHmax)>0 | as.vector(RHmin.temp)<(-50)  ))))
        if (interp_abnormal == T) {
          RHmin.temp <- ReadInput_InterpAbnormal("RHmin.temp", 
                                                 upperdata = RHmax, RHmin.temp, 
                                                 abnormal_method)
          if (is.null(RHmin.temp)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        }
      }
      
      if (length(Missing_DateIndex.daily) > 0) {
        RHmin.temp <- merge(RHmin.temp, Stdzoo, all = TRUE, 
                            fill = NA)$RHmin.temp
        if (interp_missing_days == T) {
          RHmin <- ReadInput_InterpMissing("RHmin.temp", 
                                           RHmin.temp, missing_method)
          if (is.null(RHmin)) {
            stop("More than one entry missing, please choose another interpolation method")
          }
        } else {
          RHmin <- RHmin.temp
        }
      } else {
        RHmin <- RHmin.temp
      }
    } else {
      stop("Missing variable of RHmin in input data")
    }
  }
  # Read in RH (subdaily)
  if (!"RHmax" %in% varnames & !"RHmin" %in% varnames ) {
    if ("RH" %in% varnames) {
      if (timestep != "subdaily") {
        stop("Variable RH can only be read in when date is recorded at subdaily timescale")
      }
      if ("RH" %in% (colnames(climatedata))) {
        RH.temp <- zoo(as.vector(climatedata$RH), 
                       dateagg)
        
        if ("TRUE" %in% (is.na(climatedata$RH))) {
          message("Warning: missing values in 'RH' (sub-daily relative humidity)")
          message(paste("Number of missing values in RH: ", 
                        sum(is.na(climatedata$RH))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$RH))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$RH))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$RH)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
            stop("missing data of RH exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in RH exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            RH.temp <- ReadInput_InterpMissing("RH.temp", 
                                               RH.temp, missing_method)
            if (is.null(RH.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(RH.temp) > 100|as.vector(RH.temp) <0  )) > 
            0) {
          message(paste("Number of data entries where RH has errors (RH > 100% or < 0%): ", 
                        length(which(as.vector(RH.temp) > 100|as.vector(RH.temp) < 0))))
          if (interp_abnormal == T) {
            RH.temp <- ReadInput_InterpAbnormal("RH.temp", 
                                                upperdata = NULL, RH.temp, abnormal_method)
            if (is.null(RH.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        RHmin.temp <- aggregate(RH.temp, as.Date(dateagg, 
                                                 "%d/%m/%y"), FUN = min)
        if (length(Missing_DateIndex.daily) > 0) {
          RHmin.temp <- merge(RHmin.temp, Stdzoo, all = TRUE, 
                              fill = NA)$RHmin.temp
          if (interp_missing_days == T) {
            RHmin <- ReadInput_InterpMissing("RHmin.temp", 
                                             RHmin.temp,  missing_method)
            if (is.null(RHmin)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            RHmin <- RHmin.temp
          }
        } else {
          RHmin <- RHmin.temp
        }
        
        RHmax.temp <- aggregate(RH.temp, as.Date(dateagg, 
                                                 "%d/%m/%y"), FUN = max)
        if (length(Missing_DateIndex.daily) > 0) {
          RHmax.temp <- merge(RHmax.temp, Stdzoo, all = TRUE, 
                              fill = NA)$RHmax.temp
          if (interp_missing_days == T) {
            RHmax <- ReadInput_InterpMissing("RHmax.temp", 
                                             RHmax.temp,  missing_method)
            if (is.null(RHmax)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            RHmax <- RHmax.temp
          }
        } else {
          RHmax <- RHmax.temp
        }
      } else {
        stop("Missing variable of RH in input data")
      }
          
    }   
    }
  
  # Read in u2 
  if ("u2" %in% varnames) {
    if (timestep == "daily") {
      if ("u2" %in% (colnames(climatedata))) {
        u2.temp <- zoo(as.vector(climatedata$u2), 
                          dateagg)
        if ("TRUE" %in% (is.na(climatedata$u2))) {
          message("Warning: missing values in 'u2' (daily wind speed at 2m)")
          message(paste("Number of missing values in u2: ", 
                        sum(is.na(climatedata$u2))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$u2))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$u2))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$u2)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of u2 exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in u2 exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            u2.temp <- ReadInput_InterpMissing("u2.temp", 
                                                  u2.temp, missing_method)
            if (is.null(u2.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(u2.temp) < 0  )) > 
            0) {
          message(paste("Number of day increments when u2 has errors (u2 < 0): ", 
                        length(which(as.vector(u2.temp) < 0))))
          if (interp_abnormal == T) {
            u2.temp <- ReadInput_InterpAbnormal("u2.temp", 
                                                   upperdata = NULL, u2.temp, abnormal_method)
            if (is.null(u2.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          u2.temp <- merge(u2.temp, Stdzoo, all = TRUE, 
                              fill = NA)$u2.temp
          if (interp_missing_days == T) {
            u2 <- ReadInput_InterpMissing("u2.temp", 
                                             u2.temp, missing_method)
            if (is.null(u2)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            u2 <- u2.temp
          }
        } else {
          u2 <- u2.temp
        }
      } else {
        stop("Missing variable of u2 in input data")
      }
    } else if (timestep == "subdaily") {
  # read in u2 (subdaily)  
      if ("u2" %in% (colnames(climatedata))) {
        u2.temp <- zoo(as.vector(climatedata$u2), 
                       dateagg)
        
        if ("TRUE" %in% (is.na(climatedata$u2))) {
          message("Warning: missing values in 'u2' (sub-daily wind speed at 2m)")
          message(paste("Number of missing values in u2: ", 
                        sum(is.na(climatedata$u2))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$u2))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$u2))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$u2)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
            stop("missing data of u2 exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in u2 exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            u2.temp <- ReadInput_InterpMissing("u2.temp", 
                                               u2.temp, missing_method)
            if (is.null(u2.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(u2.temp) <0  )) > 
            0) {
          message(paste("Number of data entries where u2 has errors (u2 < 0): ", 
                        length(which(as.vector(u2.temp) < 0))))
          if (interp_abnormal == T) {
            u2.temp <- ReadInput_InterpAbnormal("u2.temp", 
                                                upperdata = NULL, u2.temp, abnormal_method)
            if (is.null(u2.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        u2.temp <- aggregate(u2.temp, as.Date(dateagg, 
                                                 "%d/%m/%y"), FUN = mean)
        u2 <- u2.temp
      } else {
        stop("Missing variable of u2 in input data")
      }
    }
    
  }
  
  # Read in uz
  if ("uz" %in% varnames) {
    if (timestep == "daily") {
      if ("uz" %in% (colnames(climatedata))) {
        uz.temp <- zoo(as.vector(climatedata$uz), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$uz))) {
          message("Warning: missing values in 'uz' (daily wind speed)")
          message(paste("Number of missing values in uz: ", 
                        sum(is.na(climatedata$uz))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$uz))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$uz))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$uz)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of uz exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in uz exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            uz.temp <- ReadInput_InterpMissing("uz.temp", 
                                               uz.temp, missing_method)
            if (is.null(uz.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(uz.temp) < 0  )) > 
            0) {
          message(paste("Number of day increments when uz has errors (uz < 0): ", 
                        length(which(as.vector(uz.temp) < 0))))
          if (interp_abnormal == T) {
            uz.temp <- ReadInput_InterpAbnormal("uz.temp", 
                                                upperdata = NULL, uz.temp, abnormal_method)
            if (is.null(uz.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          uz.temp <- merge(uz.temp, Stdzoo, all = TRUE, 
                           fill = NA)$uz.temp
          if (interp_missing_days == T) {
            uz <- ReadInput_InterpMissing("uz.temp", 
                                          uz.temp, missing_method)
            if (is.null(uz)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            uz <- uz.temp
          }
        } else {
          uz <- uz.temp
        }
      } else {
        stop("Missing variable of uz in input data")
      }
    } else if (timestep == "subdaily") {
      # read in uz (subdaily)  
      if ("uz" %in% (colnames(climatedata))) {
        uz.temp <- zoo(as.vector(climatedata$uz), 
                       dateagg)
        
        if ("TRUE" %in% (is.na(climatedata$uz))) {
          message("Warning: missing values in 'uz' (sub-daily wind speed)")
          message(paste("Number of missing values in uz: ", 
                        sum(is.na(climatedata$uz))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$uz))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$uz))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$uz)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
            stop("missing data of uz exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in uz exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            uz.temp <- ReadInput_InterpMissing("uz.temp", 
                                               uz.temp, missing_method)
            if (is.null(uz.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(uz.temp) <0  )) > 
            0) {
          message(paste("Number of data entries where uz has errors (uz < 0): ", 
                        length(which(as.vector(uz.temp) < 0))))
          if (interp_abnormal == T) {
            uz.temp <- ReadInput_InterpAbnormal("uz.temp", 
                                                upperdata = NULL, uz.temp, abnormal_method)
            if (is.null(uz.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        uz.temp <- aggregate(uz.temp, as.Date(dateagg, 
                                              "%d/%m/%y"), FUN = mean)
        uz <- uz.temp
        } else {
        stop("Missing variable of uz in input data")
      }
    }
    
  }
  
  # Read in Rs
  if ("Rs" %in% varnames) {
    if (timestep == "daily") {
      if ("Rs" %in% (colnames(climatedata))) {
        Rs.temp <- zoo(as.vector(climatedata$Rs), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$Rs))) {
          message("Warning: missing values in 'Rs' (daily incoming solar radiation)")
          message(paste("Number of missing values in Rs: ", 
                        sum(is.na(climatedata$Rs))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Rs))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Rs))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Rs)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Rs exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Rs exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Rs.temp <- ReadInput_InterpMissing("Rs.temp", 
                                               Rs.temp, missing_method)
            if (is.null(Rs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(Rs.temp) < 0  )) > 
            0) {
          message(paste("Number of day increments when Rs has errors (Rs < 0): ", 
                        length(which(as.vector(Rs.temp) < 0))))
          if (interp_abnormal == T) {
            Rs.temp <- ReadInput_InterpAbnormal("Rs.temp", 
                                                upperdata = NULL, Rs.temp, abnormal_method)
            if (is.null(Rs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          Rs.temp <- merge(Rs.temp, Stdzoo, all = TRUE, 
                           fill = NA)$Rs.temp
          if (interp_missing_days == T) {
            Rs <- ReadInput_InterpMissing("Rs.temp", 
                                          Rs.temp, missing_method)
            if (is.null(Rs)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Rs <- Rs.temp
          }
        } else {
          Rs <- Rs.temp
        }
      } else {
        stop("Missing variable of Rs in input data")
      }
    } else if (timestep == "subdaily") {
      # read in Rs (subdaily)  
      if ("Rs" %in% (colnames(climatedata))) {
        Rs.temp <- zoo(as.vector(climatedata$Rs), 
                       dateagg)
        
        if ("TRUE" %in% (is.na(climatedata$Rs))) {
          message("Warning: missing values in 'Rs' (sub-daily incoming solar radiation)")
          message(paste("Number of missing values in Rs: ", 
                        sum(is.na(climatedata$Rs))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Rs))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Rs))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Rs)) >= 
              stopmissing[2]/100 * nrow(climatedata)) {
            stop("missing data of Rs exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Rs exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Rs.temp <- ReadInput_InterpMissing("Rs.temp", 
                                               Rs.temp, missing_method)
            if (is.null(Rs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(Rs.temp) <0  )) > 
            0) {
          message(paste("Number of data entries where Rs has errors (Rs < 0): ", 
                        length(which(as.vector(Rs.temp) < 0))))
          if (interp_abnormal == T) {
            Rs.temp <- ReadInput_InterpAbnormal("Rs.temp", 
                                                upperdata = NULL, Rs.temp, abnormal_method)
            if (is.null(Rs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        Rs.temp <- aggregate(Rs.temp, as.Date(dateagg, 
                                              "%d/%m/%y"), FUN = mean)
        Rs <- Rs.temp
      } else {
        stop("Missing variable of Rs in input data")
      }
    }
    
  }
  
  # Read in n
  if ("n" %in% varnames) {
    if (timestep == "daily") {
      if ("n" %in% (colnames(climatedata))) {
        n.temp <- zoo(as.vector(climatedata$n), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$n))) {
          message("Warning: missing values in 'n' (daily sunshine hours)")
          message(paste("Number of missing values in n: ", 
                        sum(is.na(climatedata$n))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$n))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$n))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$n)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of n exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in n exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            n.temp <- ReadInput_InterpMissing("n.temp", 
                                               n.temp, missing_method)
            if (is.null(n.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(n.temp) < 0  | as.vector(n.temp) > 24)) > 
            0) {
          message(paste("Number of day increments when n has errors (n < 0 or > 24 hours): ", 
                        length(which(as.vector(n.temp) < 0 | as.vector(n.temp) > 24))))
          if (interp_abnormal == T) {
            n.temp <- ReadInput_InterpAbnormal("n.temp", 
                                                upperdata = NULL, n.temp, abnormal_method)
            if (is.null(n.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          n.temp <- merge(n.temp, Stdzoo, all = TRUE, 
                           fill = NA)$n.temp
          if (interp_missing_days == T) {
            n <- ReadInput_InterpMissing("n.temp", 
                                          n.temp, missing_method)
            if (is.null(n)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            n <- n.temp
          }
        } else {
          n <- n.temp
        }
      } else {
        stop("Missing variable of n in input data")
      }
    } else if (timestep == "subdaily") {
      if ("n" %in% (colnames(climatedata))) {
        n.temp <- zoo(as.vector(climatedata$n), 
                      dateagg)
        if ("TRUE" %in% (is.na(climatedata$n))) {
          message("Warning: missing values in 'n' (daily sunshine hours)")
          message(paste("Number of missing values in n: ", 
                        sum(is.na(climatedata$n))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$n))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$n))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$n)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of n exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in n exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            n.temp <- ReadInput_InterpMissing("n.temp", 
                                              n.temp, missing_method)
            if (is.null(n.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(n.temp) < 0  | as.vector(n.temp) > 24)) > 
            0) {
          message(paste("Number of day increments when n has errors (n < 0 or > 24 hours): ", 
                        length(which(as.vector(n.temp) < 0 | as.vector(n.temp) > 24))))
          if (interp_abnormal == T) {
            n.temp <- ReadInput_InterpAbnormal("n.temp", 
                                               upperdata = NULL, n.temp, abnormal_method)
            if (is.null(n.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        n.temp <- aggregate(n.temp, as.Date(dateagg, 
                                              "%d/%m/%y"), FUN = mean)
        
        if (length(Missing_DateIndex.daily) > 0) {
          n.temp <- merge(n.temp, Stdzoo, all = TRUE, 
                          fill = NA)$n.temp
          if (interp_missing_days == T) {
            n <- ReadInput_InterpMissing("n.temp", 
                                         n.temp, missing_method)
            if (is.null(n)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            n <- n.temp
          }
        } else {
          n <- n.temp
        }
      } else {
        stop("Missing variable of n in input data")
      }
    } 
    
  }
  
  # Read in Cd
  if ("Cd" %in% varnames) {
    if (timestep == "daily") {
      if ("Cd" %in% (colnames(climatedata))) {
        C0.temp <- zoo(as.vector(climatedata$Cd), dateagg)
        if ("TRUE" %in% (is.na(climatedata$Cd))) {
          message("Warning: missing values in 'Cd' (daily cloud cover)")
          message(paste("Number of missing values in Cd: ", 
                        sum(is.na(climatedata$Cd))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Cd))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Cd))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Cd)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Cd exceeds ", stopmissing[2], 
                 "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Cd exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            C0.temp <- ReadInput_InterpMissing("C0.temp", 
                                               C0.temp, missing_method)
            if (is.null(C0.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(C0.temp) < 0)) > 0) {
          message(paste("Number of day increments when Cd has errors (Cd < 0): ", 
                        length(which(as.vector(C0.temp) < 0))))
          if (interp_abnormal == T) {
            C0.temp <- ReadInput_InterpAbnormal("C0.temp", 
                                                upperdata = NULL, C0.temp, abnormal_method)
            if (is.null(C0.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          C0.temp <- merge(C0.temp, Stdzoo, all = TRUE, fill = NA)$C0.temp
          if (interp_missing_days == T) {
            C0 <- ReadInput_InterpMissing("C0.temp",C0.temp, 
                                          missing_method)
            if (is.null(C0)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            C0 <- C0.temp
          }
        } else {
          C0 <- C0.temp
        }
        n <- constants$a_0 + constants$b_0 * C0 + constants$c_0 * 
          C0^2 + constants$d_0 * C0^3
        Cd <- C0
      } else {
        stop("Missing variable of Cd in input data")
      }
    } else if (timestep == "subdaily")  {
      if ("Cd" %in% (colnames(climatedata))) {
        C0.temp <- zoo(as.vector(climatedata$Cd), dateagg)
        if ("TRUE" %in% (is.na(climatedata$Cd))) {
          message("Warning: missing values in 'Cd' (daily cloud cover)")
          message(paste("Number of missing values in Cd: ", 
                        sum(is.na(climatedata$Cd))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Cd))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Cd))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Cd)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Cd exceeds ", stopmissing[2], 
                 "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Cd exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            C0.temp <- ReadInput_InterpMissing("C0.temp", 
                                               C0.temp, missing_method)
            if (is.null(C0.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(C0.temp) < 0)) > 0) {
          message(paste("Number of day increments when Cd has errors (Cd < 0): ", 
                        length(which(as.vector(C0.temp) < 0))))
          if (interp_abnormal == T) {
            C0.temp <- ReadInput_InterpAbnormal("C0.temp", 
                                                upperdata = NULL, C0.temp, abnormal_method)
            if (is.null(C0.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        C0.temp <- aggregate(C0.temp, as.Date(dateagg, 
                                            "%d/%m/%y"), FUN = mean)
        if (length(Missing_DateIndex.daily) > 0) {
          C0.temp <- merge(C0.temp, Stdzoo, all = TRUE, fill = NA)$C0.temp
          if (interp_missing_days == T) {
            C0 <- ReadInput_InterpMissing("C0.temp",C0.temp, 
                                          missing_method)
            if (is.null(C0)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            C0 <- C0.temp
          }
        } else {
          C0 <- C0.temp
        }
        n <- constants$a_0 + constants$b_0 * C0 + constants$c_0 * 
          C0^2 + constants$d_0 * C0^3
        Cd <- C0
      } else {
        stop("Missing variable of Cd in input data")
      }    
      }
  }
  
  # Read in Precip
  if ("Precip" %in% varnames) {
    if (timestep == "daily") {
      if ("Precip" %in% (colnames(climatedata))) {
        P.temp <- zoo(as.vector(climatedata$Precip), dateagg)
        if ("TRUE" %in% (is.na(climatedata$Precip))) {
          message("Warning: missing values in 'Precip' (daily precipitation)")
          message(paste("Number of missing values in Precip: ", 
                        sum(is.na(climatedata$Precip))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Precip))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Precip)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Precip exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Precip exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            P.temp <- ReadInput_InterpMissing("P.temp", 
                                              P.temp, missing_method)
            if (is.null(P.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(P.temp) < 0)) > 0) {
          message(paste("Number of day increments when P has errors (P < 0): ", 
                        length(which(as.vector(P.temp) < 0))))
          if (interp_abnormal == T) {
            P.temp <- ReadInput_InterpAbnormal("P.temp", 
                                               upperdata = NULL, P.temp, abnormal_method)
            if (is.null(P.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          P.temp <- merge(P.temp, Stdzoo, all = TRUE, fill = NA)$P.temp
          if (interp_missing_days == T) {
            Precip <- ReadInput_InterpMissing("P.temp", 
                                              P.temp, missing_method)
          } else {
            Precip <- P.temp
          }
        } else {
          Precip <- P.temp
        }
      } else {
        stop("Missing variable of Precip in input data")
      }
    } else if (timestep == "subdaily") {
      if ("Precip" %in% (colnames(climatedata))) {
        P.temp <- zoo(as.vector(climatedata$Precip), dateagg)
        if ("TRUE" %in% (is.na(climatedata$Precip))) {
          message("Warning: missing values in 'Precip' (sub-daily precipitation)")
          message(paste("Number of missing values in Precip: ", 
                        sum(is.na(climatedata$Precip))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Precip))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Precip)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Precip exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Precip exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            P.temp <- ReadInput_InterpMissing("P.temp", 
                                              P.temp, missing_method)
            if (is.null(P.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(P.temp) < 0)) > 0) {
          message(paste("Number of day increments when P has errors (P < 0): ", 
                        length(which(as.vector(P.temp) < 0))))
          if (interp_abnormal == T) {
            P.temp <- ReadInput_InterpAbnormal("P.temp", 
                                               upperdata = NULL, P.temp, abnormal_method)
            if (is.null(P.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        P.temp <- aggregate(P.temp, as.Date(dateagg, "%d/%m/%y"), FUN = sum)
        if (length(Missing_DateIndex.daily) > 0) {
          P.temp <- merge(P.temp, Stdzoo, all = TRUE, fill = NA)$P.temp
          if (interp_missing_days == T) {
            Precip <- ReadInput_InterpMissing("P.temp", 
                                              P.temp, missing_method)
          } else {
            Precip = P.temp
          }
        } else {
          Precip = P.temp
        }
      } else {
        stop("Missing variable of Precip in input data")
      }
      
    }
    P.monthly <- aggregate(Precip, as.yearmon, sum)
  }

  # Read in Epan
  if ("Epan" %in% varnames) {
    if (timestep == "daily") {
      if ("Epan" %in% (colnames(climatedata))) {
        Epan.temp <- zoo(as.vector(climatedata$Epan), 
                      dateagg)
        if ("TRUE" %in% (is.na(climatedata$Epan))) {
          message("Warning: missing values in 'Epan' (daily pan evaporation)")
          message(paste("Number of missing values in Epan: ", 
                        sum(is.na(climatedata$Epan))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Epan))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Epan))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Epan)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Epan exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Epan exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Epan.temp <- ReadInput_InterpMissing("Epan.temp", 
                                              Epan.temp, missing_method)
            if (is.null(Epan.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(Epan.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when Epan has errors (Epan < 0): ", 
                        length(which(as.vector(Epan.temp) < 0 ))))
          if (interp_abnormal == T) {
            Epan.temp <- ReadInput_InterpAbnormal("Epan.temp", 
                                               upperdata = NULL, Epan.temp, abnormal_method)
            if (is.null(Epan.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          Epan.temp <- merge(Epan.temp, Stdzoo, all = TRUE, 
                          fill = NA)$Epan.temp
          if (interp_missing_days == T) {
            Epan <- ReadInput_InterpMissing("Epan.temp", 
                                         Epan.temp, missing_method)
            if (is.null(Epan)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Epan <- Epan.temp
          }
        } else {
          Epan <- Epan.temp
        }
      } else {
        stop("Missing variable of Epan in input data")
      }
    } else if (timestep == "subdaily") {
      if ("Epan" %in% (colnames(climatedata))) {
        Epan.temp <- zoo(as.vector(climatedata$Epan), 
                         dateagg)
        if ("TRUE" %in% (is.na(climatedata$Epan))) {
          message("Warning: missing values in 'Epan' (daily pan evaporation)")
          message(paste("Number of missing values in Epan: ", 
                        sum(is.na(climatedata$Epan))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$Epan))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$Epan))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$Epan)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of Epan exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in Epan exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            Epan.temp <- ReadInput_InterpMissing("Epan.temp", 
                                                 Epan.temp, missing_method)
            if (is.null(Epan.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(Epan.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when Epan has errors (Epan < 0): ", 
                        length(which(as.vector(Epan.temp) < 0 ))))
          if (interp_abnormal == T) {
            Epan.temp <- ReadInput_InterpAbnormal("Epan.temp", 
                                                  upperdata = NULL, Epan.temp, abnormal_method)
            if (is.null(Epan.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        Epan.temp <- aggregate(Epan.temp, as.Date(dateagg, "%d/%m/%y"), FUN = sum)
        
        if (length(Missing_DateIndex.daily) > 0) {
          Epan.temp <- merge(Epan.temp, Stdzoo, all = TRUE, 
                             fill = NA)$Epan.temp
          if (interp_missing_days == T) {
            Epan <- ReadInput_InterpMissing("Epan.temp", 
                                            Epan.temp, missing_method)
            if (is.null(Epan)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            Epan <- Epan.temp
          }
        } else {
          Epan <- Epan.temp
        }
      } else {
        stop("Missing variable of Epan in input data")
      }    
      }
    
  }
  
  # Read in Vs
  if ("vs" %in% varnames) {
    if (timestep == "daily") {
      if ("vs" %in% (colnames(climatedata))) {
        vs.temp <- zoo(as.vector(climatedata$vs), 
                         dateagg)
        if ("TRUE" %in% (is.na(climatedata$vs))) {
          message("Warning: missing values in 'vs' (daily saturated vapour pressure)")
          message(paste("Number of missing values in vs: ", 
                        sum(is.na(climatedata$vs))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$vs))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$vs))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$vs)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of vs exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in vs exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            vs.temp <- ReadInput_InterpMissing("vs.temp", 
                                                 vs.temp, missing_method)
            if (is.null(vs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(vs.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when vs has errors (vs < 0): ", 
                        length(which(as.vector(vs.temp) < 0 ))))
          if (interp_abnormal == T) {
            vs.temp <- ReadInput_InterpAbnormal("vs.temp", 
                                                  upperdata = NULL, vs.temp, abnormal_method)
            if (is.null(vs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          vs.temp <- merge(vs.temp, Stdzoo, all = TRUE, 
                             fill = NA)$vs.temp
          if (interp_missing_days == T) {
            vs <- ReadInput_InterpMissing("vs.temp", 
                                            vs.temp, missing_method)
            if (is.null(vs)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            vs <- vs.temp
          }
        } else {
          vs <- vs.temp
        }
      } else {
        stop("Missing variable of vs in input data")
      }
    } else if (timestep == "subdaily") {
      if ("vs" %in% (colnames(climatedata))) {
        vs.temp <- zoo(as.vector(climatedata$vs), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$vs))) {
          message("Warning: missing values in 'vs' (sub-daily saturated vapour pressure)")
          message(paste("Number of missing values in vs: ", 
                        sum(is.na(climatedata$vs))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$vs))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$vs))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$vs)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of vs exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in vs exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            vs.temp <- ReadInput_InterpMissing("vs.temp", 
                                               vs.temp, missing_method)
            if (is.null(vs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(vs.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when vs has errors (vs < 0): ", 
                        length(which(as.vector(vs.temp) < 0 ))))
          if (interp_abnormal == T) {
            vs.temp <- ReadInput_InterpAbnormal("vs.temp", 
                                                upperdata = NULL, vs.temp, abnormal_method)
            if (is.null(vs.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        vs.temp <- aggregate(vs.temp, as.Date(dateagg, "%d/%m/%y"), FUN = mean)
        if (length(Missing_DateIndex.daily) > 0) {
          vs.temp <- merge(vs.temp, Stdzoo, all = TRUE, 
                           fill = NA)$vs.temp
          if (interp_missing_days == T) {
            vs <- ReadInput_InterpMissing("vs.temp", 
                                          vs.temp, missing_method)
            if (is.null(vs)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            vs <- vs.temp
          }
        } else {
          vs <- vs.temp
        }
      } else {
        stop("Missing variable of vs in input data")
      }    
      }
    
  }
  
  # Read in Va
  if ("va" %in% varnames) {
    if (timestep == "daily") {
      if ("va" %in% (colnames(climatedata))) {
        va.temp <- zoo(as.vector(climatedata$va), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$va))) {
          message("Warning: missing values in 'va' (daily average vapour pressure)")
          message(paste("Number of missing values in va: ", 
                        sum(is.na(climatedata$va))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$va))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$va))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$va)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of va exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in va exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            va.temp <- ReadInput_InterpMissing("va.temp", 
                                               va.temp, missing_method)
            if (is.null(va.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(va.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when va has errors (va < 0): ", 
                        length(which(as.vector(va.temp) < 0 ))))
          if (interp_abnormal == T) {
            va.temp <- ReadInput_InterpAbnormal("va.temp", 
                                                upperdata = NULL, va.temp, abnormal_method)
            if (is.null(va.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(Missing_DateIndex.daily) > 0) {
          va.temp <- merge(va.temp, Stdzoo, all = TRUE, 
                           fill = NA)$va.temp
          if (interp_missing_days == T) {
            va <- ReadInput_InterpMissing("va.temp", 
                                          va.temp, missing_method)
            if (is.null(va)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            va <- va.temp
          }
        } else {
          va <- va.temp
        }
      } else {
        stop("Missing variable of va in input data")
      }
    } else if (timestep == "subdaily") {
      if ("va" %in% (colnames(climatedata))) {
        va.temp <- zoo(as.vector(climatedata$va), 
                       dateagg)
        if ("TRUE" %in% (is.na(climatedata$va))) {
          message("Warning: missing values in 'va' (sub-daily average vapour pressure)")
          message(paste("Number of missing values in va: ", 
                        sum(is.na(climatedata$va))))
          message(paste("% missing data: ", signif(sum(is.na(climatedata$va))/nrow(climatedata) * 
                                                     100, digits = -3), "%"))
          x <- df <- NULL
          x <- as.numeric(!is.na(climatedata$va))
          df <- data.frame(x, zcount = NA)
          df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
          for (counter in 2:nrow(df)) {
            df$zcount[counter] <- ifelse(df$x[counter] == 
                                           0, df$zcount[counter - 1] + 1, 0)
          }
          message(paste("Maximum duration of missing data as percentage of total duration: ", 
                        signif(max(df)/nrow(climatedata) * 100, digits = -3), 
                        "%"))
          if (sum(is.na(climatedata$va)) >= stopmissing[2]/100 * 
              nrow(climatedata)) {
            stop("missing data of va exceeds ", 
                 stopmissing[2], "%, please use high quality data for calculation")
          } else {
            if (max(df)/nrow(climatedata) >= stopmissing[3]/100) {
              stop("Maximum duration of missing data in va exceeds ", 
                   stopmissing[3], "% of total data duration, please use high quality data for calculation")
            }
          }
          if (interp_missing_entries == T) {
            va.temp <- ReadInput_InterpMissing("va.temp", 
                                               va.temp, missing_method)
            if (is.null(va.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        if (length(which(as.vector(va.temp) < 0 )) > 
            0) {
          message(paste("Number of day increments when va has errors (va < 0): ", 
                        length(which(as.vector(va.temp) < 0 ))))
          if (interp_abnormal == T) {
            va.temp <- ReadInput_InterpAbnormal("va.temp", 
                                                upperdata = NULL, va.temp, abnormal_method)
            if (is.null(va.temp)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          }
        }
        va.temp <- aggregate(va.temp, as.Date(dateagg, "%d/%m/%y"), FUN = mean)
        if (length(Missing_DateIndex.daily) > 0) {
          va.temp <- merge(va.temp, Stdzoo, all = TRUE, 
                           fill = NA)$va.temp
          if (interp_missing_days == T) {
            va <- ReadInput_InterpMissing("va.temp", 
                                          va.temp, missing_method)
            if (is.null(va)) {
              stop("More than one entry missing, please choose another interpolation method")
            }
          } else {
            va <- va.temp
          }
        } else {
          va <- va.temp
        }
      } else {
        stop("Missing variable of va in input data")
      }      
    }
    
  }
#########################################
  if (!exists("Tmax")) {
    Tmax <- NULL
  }
  if (!exists("Tmin")) {
    Tmin <- NULL
  }
  if (!exists("Tdew")) {
    Tdew <- NULL
  }
  if (!exists("RHmax")) {
    RHmax <- NULL
  }
  if (!exists("RHmin")) {
    RHmin <- NULL
  }
  if (!exists("uz")) {
    uz <- NULL
  }
  if (!exists("u2")) {
    u2 <- NULL
  }
  if (!exists("Rs")) {
    Rs <- NULL
  }
  if (!exists("n")) {
    n <- NULL
  }
  if (!exists("Cd")) {
    Cd <- NULL
  }
  if (!exists("Precip")) {
    Precip <- NULL
  }
  if (!exists("P.monthly")) {
    P.monthly <- NULL
  }
  if (!exists("Epan")) {
    Epan <- NULL
  }
  if (!exists("va")) {
    va <- NULL
  }
  if (!exists("vs")) {
    vs <- NULL
  }
  data <- list(Date.daily = Date.daily, Date.monthly = Date.monthly, 
               J = J, i = i, Ndays = Ndays, Tmax = Tmax, Tdew = Tdew, Tmin = Tmin,
               RHmax = RHmax, RHmin = RHmin, u2 = u2, uz = uz,
               Rs = Rs, n = n, 
               Cd = Cd, Precip = Precip, P.monthly = P.monthly,
               Epan = Epan, va = va, vs = vs)
  invisible(data)
}
  
  ####################################
  ReadInput_InterpAbnormal <- function (varname, upperdata = NULL, var, abnormal_method) 
  {
    assign(varname, var)
    tempdata <- get(varname)
    if (grepl("Tmax", varname) | grepl("Tdew", varname) | varname == 
        "temp.temp" ) {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var) > 100 | as.vector(var) < (-20) 
        return(test)
      }
    } else if ( grepl("RHmax", varname) | varname == "RH.temp") {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var) > 100 | as.vector(var) < 0 
        return(test)
      }
    } else if (grepl("u", varname) | grepl("C", varname) | grepl("Rs", 
                                                               varname) | grepl("P", varname) | 
             grepl("Epan", varname)) {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var) < 0
        return(test)
      }
    } else if (grepl("Tmin", varname)) {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var-upperdata)>0 | as.vector(var) < (-20) 
        return(test)
      }
    } else if (grepl("RHmin", varname)) {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var-upperdata)>0 | as.vector(var) < 0 
        return(test)
      }
    } else if (grepl("n", varname)) {
      limfun <- function(varname) {
        var <- get(varname)
        test <- as.vector(var) > 24 | as.vector(var) < 0 
        return(test)
      }
    }
    
    
    tempdata <- get(varname)
    if (is.null(abnormal_method) | abnormal_method == "monthly average") {
      abnormal_method = "monthly average"
      for (m in 0:11) {
        tempdata[which(as.POSIXlt(time(tempdata))$mon == m & as.vector(limfun(varname)) == 
                   T)] = mean(tempdata[which(as.POSIXlt(time(tempdata))$mon == 
                                        m & as.vector(limfun(varname)) == F)])
      }
    }
    else if (abnormal_method == "seasonal average") {
      smonth <- rbind(c(11, 0, 1), c(2:4), c(5:7), c(8:10))
      for (s in 1:4) {
        m = c(smonth[s, ])
        tempdata[which(any(as.POSIXlt(time(tempdata))$mon == m) & 
                   as.vector(limfun(varname)) == T)] = mean(tempdata[which(as.POSIXlt(time(tempdata))$mon == 
                                                                      m & as.vector(limfun(varname)) == F)])
      }
    }
    else if (abnormal_method == "DoY average") {
      for (j in 1:366) {
        tempdata[which(as.numeric(strftime(time(tempdata), format = "%j")) == 
                   j & as.vector(limfun(varname)) == T)] = mean(tempdata[which(as.numeric(strftime(time(tempdata), 
                                                                                            format = "%j")) == j & as.vector(limfun(varname)) == 
                                                                          F)])
      }
    }
    else if (abnormal_method == "neighbour average") {
      misi <- which(is.na(tempdata[1:length(tmpdata)]))
      if (1 %in% misi) {
        tempdata[1] = tempdata[2]
      }
      if (length(tmpdata) %in% misi) {
        tempdata[length(tmpdata)] = tempdata[length(tmpdata) - 
                                               1]
      }
      for (i in misi[misi != 1 & misi != length(tmpdata)]) {
        if (!(i - 1) %in% misi) {
          if (!(i + 1) %in% misi) {
            tempdata[i] = mean(c(tempdata[i + 1], tempdata[i - 
                                                             1]))
            misi <- setdiff(misi, i)
          }
          else {
            tempdata = NULL
            misi = 0
            break
          }
        }
        if (length(misi) == 0) 
          break
      }
    }
    return(tempdata)
    message("Interpolation used to fill abnormal data entries. Method: ", 
            abnormal_method)
  }
  
  #-------------------------------------------------------------------------------------
  ReadInput_InterpMissing <- function(varname,var,missing_method) {
    assign(varname,var)
    tempdata <- get(varname)
    if (is.null(missing_method) | missing_method == "monthly average") {
      missing_method = "monthly average"
      for (m in 0:11) {
        tempdata[which(as.POSIXlt(time(tempdata))$mon == m & is.na(tempdata))] = mean(tempdata[which(as.POSIXlt(time(tempdata))$mon == m & !is.na(tempdata))])
      } 
    } else if (missing_method == "seasonal average") {
      smonth <- rbind(c(11,0,1),c(2:4),c(5:7),c(8:10))
      for (s in 1:4) {
        m = c(smonth[s,])
        tempdata[which(any(as.POSIXlt(time(tempdata))$mon == m) & is.na(tempdata))] = mean(tempdata[which(as.POSIXlt(time(tempdata))$mon == m & !is.na(tempdata))])
        
      }
    } else if (missing_method == "DoY average") {
      for (j in 1:366) {
        tempdata[which(as.numeric(strftime(time(tempdata), format = "%j")) == j & is.na(tempdata))] = mean(tempdata[which(as.numeric(strftime(time(tempdata), format = "%j")) == j & !is.na(tempdata))])
      }
    } else if (missing_method == "neighbour average") {
      #if (timestep == "daily") {
      misi <- which(is.na(tempdata[1:length(tempdata)]))
      if (1%in%misi) {
        tempdata[1] = tempdata[2]
      } 
      if (length(tempdata)%in%misi) {
        tempdata[length(tempdata)] = tempdata[length(tempdata)-1]
      } 
      for (i in misi[misi!=1&misi!=length(tempdata)]) {
        if (!(i-1)%in%misi) { # means i is the start of missing value
          if (!(i+1)%in%misi) { # means i is the finish of missing value i.e. only one missing
            tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
            misi <- setdiff(misi,i)
          } else { # means i is not the finish of missing value i.e. more than one missing
            #stop("More than one entry missing, please choose another interpolation method")
            tempdata = NULL
            misi  = 0
            break
            #fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
            #nmi <- fi-i+1
            #tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
            #misi <- setdiff(misi,i:fi)
          }
        } 
        if (length(misi) == 0) break
      }
      #} else if (timestep == "subdaily") {
      # misi <- which(is.na(tempdata[1:length(Date.subdaily)]))
      #  if (1%in%misi) {
      #  tempdata[1] = tempdata[2]
      #} 
      #if (length(Date.subdaily)%in%misi) {
      #  tempdata[length(Date.subdaily)] = tempdata[length(Date.subdaily)-1]
      #} 
      #for (i in misi[misi!=1&misi!=length(Date.subdaily)]) {
      #  if (!(i-1)%in%misi) { # means i is the start of missing value
      #    if (!(i+1)%in%misi) { # means i is the finish of missing value i.e. only one missing
      #      tempdata[i] = mean(c(tempdata[i+1],tempdata[i-1]))
      #      misi <- setdiff(misi,i)
      #    } else { # means i is not the finish of missing value i.e. more than one missing
      #      fi <- c(misi[misi>i & !(misi+1)%in%misi])[1]
      #      nmi <- fi-i+1
      #      tempdata[i:fi] <- as.vector(tempdata[i-1])+(as.vector(tempdata[fi+1])-as.vector(tempdata[i-1]))/(nmi+1)*(1:nmi)
      #      misi <- setdiff(misi,i:fi)
      #    }
      #  } 
      #  if (length(misi) == 0) break
    }
    
    
    return(tempdata)
    message("Interpolation used to fill missing data entries
          . Method: ", missing_method)
  }
  
  #####################################
  
  ReadOBSEvaporation <- function(E_OBS, data) {
    
    # Load evaporation observations and convert to zoo object
    Date.OBS <- strptime(paste(E_OBS$Day, "/", E_OBS$Month, "/", E_OBS$Year, sep=""), "%d/%m/%Y")
    if (E_OBS$Day[1] == E_OBS$Day[2]) {
      Date.OBS <- as.yearmon(Date.OBS)
      E_obs <- zoo(E_OBS$EVAP.Obs, Date.OBS)
      
      # Aggregation
      E_obs.Daily <- NULL
      E_obs.Monthly <- E_obs
      E_obs.Annual <- aggregate(E_obs.Monthly, floor(as.numeric(as.yearmon(Date.OBS, "%y"))), FUN = sum)
      
      # Average
      E_obs.MonthlyAve <- E_obs.AnnualAve <- NULL
      E_obs.MonthlyAve.temp <- E_obs.Monthly/data$ndays
      for (mon in min(as.POSIXlt(Date.OBS)$mon):max(as.POSIXlt(Date.OBS)$mon)){
        i = mon - min(as.POSIXlt(Date.OBS)$mon) + 1
        E_obs.MonthlyAve[i] <- mean(E_obs.MonthlyAve.temp[as.POSIXlt(Date.OBS)$mon== mon])
      }
      for (year in min(as.POSIXlt(Date.OBS)$year):max(as.POSIXlt(Date.OBS)$year)){
        i = year - min(as.POSIXlt(Date.OBS)$year) + 1
        E_obs.AnnualAve[i] <- mean(E_obs.MonthlyAve.temp[as.POSIXlt(Date.OBS)$year== year])
      }
    } else {
      Date.OBS <- unique(as.Date(Date.OBS, "%d/%m/%y"))
      E_obs <- zoo(E_OBS$EVAP.Obs, Date.OBS)
      
      # Aggregation
      E_obs.Daily <- E_obs
      E_obs.Monthly <- aggregate(E_obs, as.yearmon(Date.OBS, "%m/%y"), FUN = sum)
      E_obs.Annual <- aggregate(E_obs.Daily, floor(as.numeric(as.yearmon(Date.OBS, "%y"))), FUN = sum)
      
      # Average
      E_obs.MonthlyAve <- E_obs.AnnualAve <- NULL
      for (mon in min(as.POSIXlt(Date.OBS)$mon):max(as.POSIXlt(Date.OBS)$mon)){
        i = mon - min(as.POSIXlt(Date.OBS)$mon) + 1
        E_obs.MonthlyAve[i] <- mean(E_obs.Daily[as.POSIXlt(Date.OBS)$mon== mon])
      }
      for (year in min(as.POSIXlt(Date.OBS)$year):max(as.POSIXlt(Date.OBS)$year)){
        i = year - min(as.POSIXlt(Date.OBS)$year) + 1
        E_obs.AnnualAve[i] <- mean(E_obs.Daily[as.POSIXlt(Date.OBS)$year== year])
      }
    }
    
    
    OBS <- list(Date.OBS=Date.OBS, E_obs.Daily=E_obs.Daily, E_obs.Monthly=E_obs.Monthly, E_obs.Annual=E_obs.Annual, E_obs.MonthlyAve=E_obs.MonthlyAve, E_obs.AnnualAve=E_obs.AnnualAve)
    return(OBS)
  }
  
  
  