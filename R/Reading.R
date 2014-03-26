# For global variable 'funname'
if(getRversion() >= "2.15.1")  utils::globalVariables(c("funname"))


# A Function for reading, checking and doing basic calculation from data for Evapotranspiration functions #
# Timestep - daily

ReadInputs <- function(climatedata, constants, stopmissing) {
  
  # Checking if all data required are available, give error message if not
  
  # Check if data 'Year', 'Month' and 'Day' exist
  if ("Year" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Year'")
  }
  if ("Month" %in% (colnames(climatedata)) == FALSE) {
    stop("missing data of 'Month'")
  }
  if ("Day" %in% (colnames(climatedata)) == FALSE) {  
    stop("missing data of 'Day'")     
  }
  
  # Date
  Date.subdaily <- strptime(paste(climatedata$Day, "/", climatedata$Month, "/", climatedata$Year, " ", climatedata$Hour, sep=""), "%d/%m/%Y %H")
  Date.daily <- unique(as.Date(Date.subdaily, "%d/%m/%y"))
  Date.monthly <- unique(as.yearmon(Date.subdaily, "%d/%m/%y"))
  
  # Julian day
  J.temp <- zoo(Date.subdaily$yday+1,as.Date(Date.subdaily)) # Julian calendar day; 1 added in to that first Jan = day 1
  J <- aggregate(J.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  
  # Number of month
  i.temp  <- unique(as.yearmon(Date.daily, "%m/%y"))
  i  <- (i.temp - trunc(i.temp))*12 + 1    #Month 
  
  # Number of days in a month
  ndays.temp <- zoo(climatedata$Day, as.Date(Date.subdaily))
  ndays.temp <- aggregate(ndays.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  ndays <- aggregate(ndays.temp, as.yearmon(Date.daily, "%m/%y"), FUN = max)
  
  # check acceptable % missing data
  if (is.na(as.numeric(stopmissing))) {
    stop("Please use a numeric value for the maximum allowable percentage of missing data")
  }
  if (!is.na(as.numeric(stopmissing))) {
    if (length(stopmissing)!=2) {
      stop("Please input a vector of length 2 for argument 'stopmissing'")
    }
  }
  if (!is.na(as.numeric(stopmissing))) {
    if (as.numeric(stopmissing) < 1 | as.numeric(stopmissing) > 99) {
      stop("Please use a value between 1 and 99 for the maximum allowable percentage of missing data")
    }
  } 
  
  message(paste("The maximum acceptable percentage of missing data is", stopmissing, "%"))
  
  # Check if data 'Tmax.daily' and 'Tmin.daily' exist
  
  if ("Tmax.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Tmax.daily))) {
      message("Warning: missing values in 'Tmax.daily' (daily maximum temperature)")
      message(paste("Number of missing values in Tmax.daily: ", sum(is.na(climatedata$Tmax.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmax.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Tmax.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Tmax.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Tmax.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Tmax.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    Tmax.temp <- zoo(climatedata$Tmax.daily, as.Date(Date.subdaily))
    Tmax <- aggregate(Tmax.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tmax has errors: ", sum(Tmax>100)))
    message("Monthly averages have been calculated to adjust data with error")
    for (m in 0:11) {
      Tmax[as.POSIXlt(time(Tmax))$mon==m & as.numeric(Tmax)>100] = mean(Tmax[as.POSIXlt(time(Tmax))$mon==m & as.numeric(Tmax)<100])
      Tmax[as.POSIXlt(time(Tmax))$mon==m & is.na(Tmax)] = mean(Tmax[as.POSIXlt(time(Tmax))$mon==m & !is.na(Tmax)])
    }
  } else {
    if ("Temp.subdaily" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'Temp.daily'(daily maximum temperature), calculated from subdaily 'Temp.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$Temp.subdaily))) {
        message("Warning: missing values in 'Temp.subdaily'")
        message(paste("Number of missing values in Temp.subdaily: ", sum(is.na(climatedata$Temp.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Temp.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
        
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$Temp.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for(i in 2:nrow(df)) {
          df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
        
        if (sum(is.na(climatedata$Temp.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
          stop("missing data of Temp.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
            stop("Maximum duration of missing data in Temp.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
          }
          message("Monthly averages have been calculated to fill missing data entries")
        }
        temp.temp <- zoo(climatedata$Temp.subdaily, as.Date(Date.subdaily))
        for (m in 0:11) {
          temp.temp[as.POSIXlt(time(temp.temp))$mon==m & is.na(temp.temp)] = mean(temp.temp[as.POSIXlt(time(temp.temp))$mon==m & !is.na(temp.temp)])
        }
        Temp <- aggregate(temp.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = max)
      }
    } else {
      Tmax <- NULL
    }
  }
  
  if ("Tmin.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Tmin.daily))) {
      message("Warning: missing values in 'Tmin.daily' (daily minimum temperature)")
      message(paste("Number of missing values in Tmin.daily: ", sum(is.na(climatedata$Tmin.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmin.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Tmin.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Tmin.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Tmin.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Tmin.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    Tmin.temp <- zoo(climatedata$Tmin.daily, as.Date(Date.subdaily))
    Tmin <- aggregate(Tmin.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tmin has errors: ", sum(Tmin>100)))
    message("Monthly averages have been calculated to adjust data with error")
    for (m in 0:11) {
      Tmin[as.POSIXlt(time(Tmin))$mon==m & as.numeric(Tmin)>100] = mean(Tmin[as.POSIXlt(time(Tmin))$mon==m & as.numeric(Tmin)<100])
      Tmin[as.POSIXlt(time(Tmin))$mon==m & is.na(Tmin)] = mean(Tmin[as.POSIXlt(time(Tmin))$mon==m & !is.na(Tmin)])
    }      
  } else if ("Temp.subdaily" %in% (colnames(climatedata))) {
    message("Warning: missing data of 'Tmin.daily'(daily minimum temperature), calculated from subdaily 'Temp.subdaily'")
    if ("TRUE" %in% (is.na(climatedata$Temp.subdaily))) {
      message("Warning: missing values in 'Temp.subdaily'")
      message(paste("Number of missing values in Temp.subdaily: ", sum(is.na(climatedata$Temp.subdaily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Temp.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Temp.subdaily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Temp.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Temp.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Temp.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
      temp.temp <- zoo(climatedata$Temp.subdaily, as.Date(Date.subdaily))
      for (m in 0:11) {
        temp.temp[as.POSIXlt(time(temp.temp))$mon==m & is.na(temp.temp)] = mean(temp.temp[as.POSIXlt(time(temp.temp))$mon==m & !is.na(temp.temp)])
      }
      Tmin <- aggregate(temp.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = min)
    }
  } else {
    Tmin <- NULL
  }
  
  Ta <- NULL
  
  # Check if data 'uz.subdaily' or 'u2.subdaily' exists
  if ("u2.subdaily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$u2.subdaily))) {
      message("Warning: missing values in 'u2.subdaily'")
      message(paste("Number of missing values in u2.subdaily: ", sum(is.na(climatedata$u2.subdaily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$u2.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$u2.subdaily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$u2.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of u2.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in u2.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    u2.temp <- zoo(climatedata$u2.subdaily*1000/3600, as.Date(Date.subdaily))
    for (m in 0:11) {
      u2.temp[as.POSIXlt(time(u2.temp))$mon==m &  as.numeric(is.na(u2.temp))] = mean(u2.temp[as.POSIXlt(time(u2.temp))$mon==m &  as.numeric(!is.na(u2.temp))]) # Changing to monthly mean (once again doesn't affect large portion of the sample)
    }
    u2 <- aggregate(u2.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  } else if ("uz.subdaily" %in% (colnames(climatedata))) {
    if ("TRUE" %in% (is.na(climatedata$uz.subdaily))) {
      message("Warning: missing data of 'u2.subdaily', calculated from 'uz.subdaily")
      u2 <- NULL
      message(paste("Number of missing values in uz.subdaily: ", sum(is.na(climatedata$uz.subdaily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$uz.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$uz.subdaily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$uz.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of uz.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in uz.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    uz.temp <- zoo(climatedata$uz.subdaily*1000/3600, as.Date(Date.subdaily)) 
    for (m in 0:11) {
      uz.temp[as.POSIXlt(time(uz.temp))$mon==m &  as.numeric(is.na(uz.temp))] = mean(uz.temp[as.POSIXlt(time(uz.temp))$mon==m &  as.numeric(!is.na(uz.temp))]) # Changing to monthly mean (once again doesn't affect large portion of the sample)
    }
    uz <- aggregate(uz.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean) 
  } else {
    u2 <- NULL
    uz <- NULL
  }
  
  # Check if data 'Rs.daily' exists
  if ('Rs.daily' %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Rs.daily))) {
      message("Warning: missing values in 'Tdew.subdaily'")
      message(paste("Number of missing values in Tdew.subdaily: ", sum(is.na(climatedata$Rs.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Rs.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Rs.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Rs.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Rs.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Rs.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    Rs.temp <- zoo(climatedata$Rs.daily, as.Date(Date.subdaily))
    Rs <- aggregate(Rs.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message("Monthly averages have been calculated to adjust data with error")
    for (m in 0:11) {
      Rs[as.POSIXlt(time(Rs))$mon==m & is.na(Rs)] = mean(Rs[as.POSIXlt(time(Rs))$mon==m & !is.na(Rs)])
    }
  } else {
    Rs <- NULL
  }
  # Check if data 'n.daily' exists
  if ("n.daily" %in% (colnames(climatedata))) {
    if ("TRUE" %in% (is.na(climatedata$n.daily))) {
      message("Warning: missing values in 'n' (daily sunshine hours)")
      message(paste("Number of missing values in n: ", sum(is.na(climatedata$n.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$n.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$n.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$n.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of n.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in n.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    n.temp <- zoo(climatedata$n.daily, as.Date(Date.subdaily))  
    for (m in 0:11) {
      n.temp[as.POSIXlt(time(n.temp))$mon==m &  as.numeric(is.na(n.temp))] = mean(n.temp[as.POSIXlt(time(n.temp))$mon==m &  as.numeric(!is.na(n.temp))]) # Changing to monthly mean 
    }
    n <- aggregate(n.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  } else {
    n <- NULL
  } 
  
  # Check if data 'Cd.daily' exists
  if ("Cd.daily" %in% (colnames(climatedata))) {
    if ("TRUE" %in% (is.na(climatedata$Cd.daily))) {
      message("Warning: missing values in 'Cd.daily' (daily cloud cover)")
      message(paste("Number of missing values in Cd.daily: ", sum(is.na(climatedata$Cd.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Cd.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Cd.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Cd.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Cd.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Cd.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    C0.temp <- zoo(climatedata$Cd.daily, as.Date(Date.subdaily)) 
    for (m in 0:11) {
      C0.temp[as.POSIXlt(time(C0.temp))$mon==m &  as.numeric(is.na(C0.temp))] = mean(C0.temp[as.POSIXlt(time(C0.temp))$mon==m &  as.numeric(!is.na(C0.temp))]) # Changing to monthly mean 
    }
    C0 <- aggregate(C0.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    n <- constants$a_0 + constants$b_0*C0 + constants$c_0*C0^2 + constants$d_0*C0^3 # calculation of sunshine hours (h) based on cloud cover (oktas) (S3.10)
  } 
  
  
  # Check if data 'Precip.daily' exists
  if ("Precip.daily" %in% (colnames(climatedata))) {
    if ("TRUE" %in% (is.na(climatedata$Precip.daily))) {
      message("Warning: missing values in 'Precip.daily' (daily precipitation)")
      message(paste("Number of missing values in Precip.daily: ", sum(is.na(climatedata$Precip.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Precip.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Precip.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Precip.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Precip.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    P.temp <- zoo(climatedata$Precip.daily, as.Date(Date.subdaily)) 
    for (m in 0:11) {
      P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(is.na(P.temp))] = mean(P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(!is.na(P.temp))]) # Changing to monthly mean 
    }
    Precip <- aggregate(P.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    P.monthly <- aggregate(Precip, as.yearmon(Date.daily, "%m/%y"),sum)
    Cd.temp <- numeric(length(P.monthly))
    for (m in 1:length(P.monthly)) {
      if (P.monthly[m] >= 1) {
        Cd.temp[m] <- 1 + 0.5*log10(P.monthly[m]) + (log10(P.monthly[m]))^2 # calculation of cloudiness (number of tenths of sky covered by cloud) based on monthly precipitation (mm) (S3.12)
      } else {
        Cd.temp[m] <- 1 # calculation of cloudiness (number of tenths of sky covered by cloud) based on monthly precipitation (mm) (S3.12)
      }
    }
    Cd.temp <- zoo(Cd.temp, as.Date(Date.daily))
    Cd <- Precip
    for (m in 1:length(Cd)) {
      Cd[as.yearmon(time(Cd)) == as.yearmon(time(Cd.temp))[m]] <- Cd.temp[m]
    }
  } else {
    Cd <- NULL
  } 
  
  
  # Check if data 'Precip.daily' exist
  if ("Precip.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Precip.daily))) {
      message("Warning: missing values in 'Precip.daily' (daily precipitation)")
      message(paste("Number of missing values in Precip.daily: ", sum(is.na(climatedata$Precip.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Precip.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Precip.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Precip.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Precip.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Precip.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    P.temp <- zoo(climatedata$Precip.daily, as.Date(Date.subdaily)) 
    for (m in 0:11) {
      P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(is.na(P.temp))] = mean(P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(!is.na(P.temp))]) # Changing to monthly mean 
    }
    Precip <- aggregate(P.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  } else {
    Precip <- NULL
  }
  
  # Check if data 'Epan.daily' exist
  if ("Epan.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Epan.daily))) {
      message("Warning: missing values in 'Epan.daily'")
      message(paste("Number of missing values in Epan.daily: ", sum(is.na(climatedata$Epan.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Epan.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Epan.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Epan.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Epan.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Epan.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    Epan.temp <- zoo(climatedata$Epan.daily, as.Date(Date.subdaily)) 
    for (m in 0:11) {
      Epan.temp[as.POSIXlt(time(Epan.temp))$mon==m &  as.numeric(is.na(Epan.temp))] = mean(Epan.temp[as.POSIXlt(time(Epan.temp))$mon==m &  as.numeric(!is.na(Epan.temp))]) # Changing to monthly mean 
    }
    Epan <- aggregate(P.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean) 
  } else {
    Epan <- NULL
  }
  
  # Check if data 'RHmax.daily' and 'RHmin.daily' exist
  if ("RHmax.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$RHmax.daily))) {
      message("Warning: missing values in 'RHmax.daily' (daily maximum temperature)")
      message(paste("Number of missing values in RHmax.daily: ", sum(is.na(climatedata$RHmax.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmax.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$RHmax.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$RHmax.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of RHmax.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in RHmax.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    RHmax.temp <- zoo(climatedata$RHmax.daily, as.Date(Date.subdaily))
    RHmax <- aggregate(RHmax.temp, as.Date(Date.subdaily, "%d/%m/%y") ,FUN = mean) 
    for (m in 0:11) {
      RHmax[as.POSIXlt(time(RHmax))$mon==m & is.na(RHmax)] = mean(RHmax[as.POSIXlt(time(RHmax))$mon==m & !is.na(RHmax)])
    }
    message(paste("Number of days increments when RHmax has errors: ", sum(RHmax>100)))
    message("'RHmax.daily' with values > 100% has been adjusted to 100%")
    RHmax[RHmax>100] = 100
  } else {
    if ("RH.subdaily" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'RHmax.daily'(daily maximum relative humidity), calculated from subdaily 'RH.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$RH.subdaily))) {
        message("Warning: missing values in 'RH.subdaily'")
        message(paste("Number of missing values in RH.subdaily: ", sum(is.na(climatedata$RH.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RH.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
        
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RH.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for(i in 2:nrow(df)) {
          df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
        
        if (sum(is.na(climatedata$RH.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
          stop("missing data of RH.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
            stop("Maximum duration of missing data in RH.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
          }
          message("Monthly averages have been calculated to fill missing data entries")
        }
      }
      RH.temp <- zoo(climatedata$RH.subdaily, as.Date(Date.subdaily))
      RHmax <- aggregate(RH.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = max)
      for (m in 0:11) {
        RHmax[as.POSIXlt(time(RHmax))$mon==m & is.na(RHmax)] = mean(RHmax[as.POSIXlt(time(RHmax))$mon==m & !is.na(RHmax)])
      }
      message(paste("Number of days increments when RHmax has errors: ", sum(RHmax>100)))
      message("'RHmax.daily' with values > 100% has been adjusted to 100%")
      RHmax[RHmax>100] = 100
    } else {
      RHmax <- NULL
    }
  }
  
  if ("RHmin.daily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$RHmin.daily))) {
      message("Warning: missing values in 'RHmin.daily' (daily minimum temperature)")
      message(paste("Number of missing values in RHmin.daily: ", sum(is.na(climatedata$RHmin.daily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmin.daily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$RHmin.daily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$RHmin.daily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of RHmin.daily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in RHmin.daily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
    RHmin.temp <- zoo(climatedata$RHmin.daily, as.Date(Date.subdaily))
    RHmin <- aggregate(RHmin.temp, as.Date(Date.subdaily, "%d/%m/%y") ,FUN = mean) 
    for (m in 0:11) {
      RHmin[as.POSIXlt(time(RHmin))$mon==m & is.na(RHmin)] = mean(RHmin[as.POSIXlt(time(RHmin))$mon==m & !is.na(RHmin)])
    }
    message(paste("Number of days increments when RHmin has errors: ", sum(RHmin>RHmax)))
    message("'RHmin.daily' with values > 'RHmax.daily' has been adjusted to '0.9*RHmax'")
    RHmin[RHmin>=RHmax] = 0.9*RHmax[RHmin>=RHmax] # setting upper bound on minimum RH to ensure that ea < es for all es,ea. 
  } else {
    if ("RH.subdaily" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'RHmin.daily'(daily minimum relative humidity), calculated from subdaily 'RH.subdaily'")
      if ("TRUE" %in% (is.na(climatedata$RH.subdaily))) {
        message("Warning: missing values in 'RH.subdaily'")
        message(paste("Number of missing values in RH.subdaily: ", sum(is.na(climatedata$RH.subdaily))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$RH.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
        x <- df <- NULL
        x <- as.numeric(!is.na(climatedata$RH.subdaily))
        df <- data.frame(x, zcount = NA)
        df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
        for(i in 2:nrow(df)) {
          df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
        }
        message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
        
        if (sum(is.na(climatedata$RH.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
          stop("missing data of RH.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
        } else {
          if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
            stop("Maximum duration of missing data in RH.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
          }
          message("Monthly averages have been calculated to fill missing data entries")
        }
      }
      RH.temp <- zoo(climatedata$RH.subdaily, as.Date(Date.subdaily))
      RHmin <- aggregate(RH.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = min)
      for (m in 0:11) {
        RHmin[as.POSIXlt(time(RHmin))$mon==m & is.na(RHmin)] = mean(RHmin[as.POSIXlt(time(RHmin))$mon==m & !is.na(RHmin)])
      }
      message(paste("Number of days increments when RHmin has errors: ", sum(RHmin>RHmax)))
      message("'RHmax' with values > 'RHmax' has been adjusted to '0.9*RHmax'")
      RHmin[RHmin>=RHmax] = 0.9*RHmax[RHmin>=RHmax] # setting upper bound on minimum RH to ensure that ea < es for all es,ea. 
    } else {
      RHmin <- NULL
    }
  }
  
  # Check if data 'Tdew.subdaily' exists
  if ("Tdew.subdaily" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Tdew.subdaily))) {
      message("Warning: missing values in 'Tdew.subdaily'")
      message(paste("Number of missing values in Tdew.subdaily: ", sum(is.na(climatedata$Tdew.subdaily))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tdew.subdaily))/nrow(climatedata) * 100, digits = -3), "%"))
      x <- df <- NULL
      x <- as.numeric(!is.na(climatedata$Tdew.subdaily))
      df <- data.frame(x, zcount = NA)
      df$zcount[1] <- ifelse(df$x[1] == 0, 1, 0)
      for(i in 2:nrow(df)) {
        df$zcount[i] <- ifelse(df$x[i] == 0, df$zcount[i - 1] + 1, 0)
      }
      message(paste("Maximum duration of missing data as percentage of total duration: ", signif(max(df)/nrow(climatedata) * 100, digits = -3), "%"))
      
      if (sum(is.na(climatedata$Tdew.subdaily)) >= stopmissing[1]/100 * nrow(climatedata)) {
        stop("missing data of Tdew.subdaily exceeds ", stopmissing[1], "%, please use high quality data for calculation")
      } else {
        if (max(df)/nrow(climatedata) >= stopmissing[2]/100) {
          stop("Maximum duration of missing data in Tdew.subdaily exceeds ", stopmissing[2], "% of total data duration, please use high quality data for calculation")
        }
        message("Monthly averages have been calculated to fill missing data entries")
      }
    }
      
    Tdew.temp <- zoo(climatedata$Tdew.subdaily, as.Date(Date.subdaily))
    Tdew <- aggregate(Tdew.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tdew has errors: ", sum(Tdew>100)))
    message("Monthly averages have been calculated to adjust data with error")
    for (m in 0:11) {
      Tdew[as.POSIXlt(time(Tdew))$mon==m & is.na(Tdew)] = mean(Tdew[as.POSIXlt(time(Tdew))$mon==m & !is.na(Tdew)])
      Tdew[as.POSIXlt(time(Tdew))$mon==m & as.numeric(Tdew)>100] = mean(Tdew[as.POSIXlt(time(Tdew))$mon==m & as.numeric(Tdew)<100])
    }
  } else {
    if ("Vp.subdaily" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'Tdew.subdaily', calculated from 'Vp.subdaily'")
      Tdew <- NULL
      if ("TRUE" %in% (is.na(climatedata$Vp.subdaily))) {
        message("Warning: missing values in 'Vp.subdaily'")
      }
    } else {
      Tdew <- NULL
    }
  }
  
  #-------------------------------------------------------------------------------------
  
  data <- list(Date.daily=Date.daily, Date.monthly=Date.monthly, J=J, i=i, ndays=ndays, Tmax=Tmax, Tmin=Tmin, u2=u2, uz=uz, 
               Rs=Rs,n=n, Cd=Cd, Precip=Precip, Epan=Epan, RHmax=RHmax, RHmin=RHmin, Tdew=Tdew)
  return(data)
}

#-------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------