# For global variable 'funname'
if(getRversion() >= "2.15.1")  utils::globalVariables(c("funname"))


# A Function for reading, checking and doing basic calculation from data for Evapotranspiration functions #
# Timestep - daily

 Reading <- function(climatedata, constants, stopmissing) {
  
  # Checking if all data required is available, give error message if not
  
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
      if (as.numeric(stopmissing) < 1 | as.numeric(stopmissing) > 99) {
        stop("Please use a value between 1 and 99 for the maximum allowable percentage of missing data")
      }
    } 

  message(paste("The maximum acceptable percentage of missing data is", stopmissing, "%"))
  
  # Check if data 'Tmax' and 'Tmin' exist

  if ("Tmax" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Tmax))) {
      message("Warning: missing values in 'Tmax' (daily maximum temperature)")
      message(paste("Number of missing values in Tmax: ", sum(is.na(climatedata$Tmax))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmax))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Tmax)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Tmax exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    Tmax.temp <- zoo(climatedata$Tmax, Date.subdaily)  
    Tmax <- aggregate(Tmax.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tmax has errors: ", sum(Tmax>100)))
    message("Interpolation is performed to replace data with error, users should be aware of the associated risk")
    for (m in 0:11) {
      Tmax[as.POSIXlt(time(Tmax))$mon==m & as.numeric(Tmax)>100] = mean(Tmax[as.POSIXlt(time(Tmax))$mon==m & as.numeric(Tmax)<100])
      Tmax[as.POSIXlt(time(Tmax))$mon==m & is.na(Tmax)] = mean(Tmax[as.POSIXlt(time(Tmax))$mon==m & !is.na(Tmax)])
    }
  } else {
    if ("Air.temperature.in.Degrees.C" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'Tmax'(daily maximum temperature), calculated from subdaily 'Air.temperature.in.Degrees.C'")
      if ("TRUE" %in% (is.na(climatedata$Air.temperature.in.Degrees.C))) {
        message("Warning: missing values in 'Air.temperature.in.Degrees.C'")
        message(paste("Number of missing values in Air.temperature.in.Degrees.C: ", sum(is.na(climatedata$Air.temperature.in.Degrees.C))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Air.temperature.in.Degrees.C))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$Air.temperature.in.Degrees.C)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of Air.temperature.in.Degrees.C exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
        temp.temp <- zoo(climatedata$Air.temperature.in.Degrees.C, Date.subdaily)
        for (m in 0:11) {
          temp.temp[as.POSIXlt(time(temp.temp))$mon==m & is.na(temp.temp)] = mean(temp.temp[as.POSIXlt(time(temp.temp))$mon==m & !is.na(temp.temp)])
        }
        Tmax <- aggregate(temp.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = max)
    }
      } else {
      Tmax <- NULL
    }
  }
  
  if ("Tmin" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Tmin))) {
      message("Warning: missing values in 'Tmin' (daily minimum temperature)")
      message(paste("Number of missing values in Tmin: ", sum(is.na(climatedata$Tmin))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Tmin))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Tmin)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Tmin exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    Tmin.temp <- zoo(climatedata$Tmin, Date.subdaily)
    Tmin <- aggregate(Tmin.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tmin has errors: ", sum(Tmin>100)))
    message("Interpolation is performed to replace data with error, users should be aware of the associated risk")
    for (m in 0:11) {
      Tmin[as.POSIXlt(time(Tmin))$mon==m & as.numeric(Tmin)>100] = mean(Tmin[as.POSIXlt(time(Tmin))$mon==m & as.numeric(Tmin)<100])
      Tmin[as.POSIXlt(time(Tmin))$mon==m & is.na(Tmin)] = mean(Tmin[as.POSIXlt(time(Tmin))$mon==m & !is.na(Tmin)])
    }      
  } else if ("Air.temperature.in.Degrees.C" %in% (colnames(climatedata))) {
      message("Warning: missing data of 'Tmin'(daily minimum temperature), calculated from subdaily 'Air.temperature.in.Degrees.C'")
      if ("TRUE" %in% (is.na(climatedata$Air.temperature.in.Degrees.C))) {
        message("Warning: missing values in 'Air.temperature.in.Degrees.C'")
        message(paste("Number of missing values in Air.temperature.in.Degrees.C: ", sum(is.na(climatedata$Air.temperature.in.Degrees.C))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Air.temperature.in.Degrees.C))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$Air.temperature.in.Degrees.C)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of Air.temperature.in.Degrees.C exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
        temp.temp <- zoo(climatedata$Air.temperature.in.Degrees.C, Date.subdaily)
        for (m in 0:11) {
          temp.temp[as.POSIXlt(time(temp.temp))$mon==m & is.na(temp.temp)] = mean(temp.temp[as.POSIXlt(time(temp.temp))$mon==m & !is.na(temp.temp)])
        }
        Tmin <- aggregate(temp.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = min)
    }
      } else {
      Tmin <- NULL
    }
  
  Ta <- NULL
  
  # Check if data 'Wind.speed.measured.in.km.h' exists
  if ("Wind.speed.at.2m.measured.in.km.h" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Wind.speed.at.2m.measured.in.km.h))) {
      message("Warning: missing values in 'Wind.speed.at.2m.measured.in.km.h'")
      message(paste("Number of missing values in Wind.speed.at.2m.measured.in.km.h: ", sum(is.na(climatedata$Wind.speed.at.2m.measured.in.km.h))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Wind.speed.at.2m.measured.in.km.h))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Wind.speed.at.2m.measured.in.km.h)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Wind.speed.at.2m.measured.in.km.h exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    u2.temp <- zoo(climatedata$Wind.speed.at.2m.measured.in.km.h*1000/3600, Date.subdaily) 
    for (m in 0:11) {
      u2.temp[as.POSIXlt(time(u2.temp))$mon==m &  as.numeric(is.na(u2.temp))] = mean(u2.temp[as.POSIXlt(time(u2.temp))$mon==m &  as.numeric(!is.na(u2.temp))]) # Changing to monthly mean (once again doesn't affect large portion of the sample)
    }
    u2 <- aggregate(u2.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    } else if ("Wind.speed.measured.in.km.h" %in% (colnames(climatedata))) {
      if ("TRUE" %in% (is.na(climatedata$Wind.speed.measured.in.km.h))) {
      message("Warning: missing data of 'Wind.speed.measured.in.km.h', calculated from 'Wind.speed.measured.in.km.h")
      u2 <- NULL
      message(paste("Number of missing values in Wind.speed.measured.in.km.h: ", sum(is.na(climatedata$Wind.speed.measured.in.km.h))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Wind.speed.measured.in.km.h))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Wind.speed.measured.in.km.h)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Wind.speed.measured.in.km.h exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      uz.temp <- zoo(climatedata$Wind.speed.measured.in.km.h*1000/3600, Date.subdaily)  
      for (m in 0:11) {
        uz.temp[as.POSIXlt(time(uz.temp))$mon==m &  as.numeric(is.na(uz.temp))] = mean(uz.temp[as.POSIXlt(time(uz.temp))$mon==m &  as.numeric(!is.na(uz.temp))]) # Changing to monthly mean (once again doesn't affect large portion of the sample)
      }
      uz <- aggregate(uz.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean) 
      } else {
      u2 <- NULL
      uz <- NULL
    }
  
  # Check if data 'Solar.radiation.in.MJ.per.sqm.per.day' exists
  if ('Solar.radiation.in.MJ.per.sqm.per.day' %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Solar.radiation.in.MJ.per.sqm.per.day))) {
      message("Warning: missing values in 'Dew.point.temperature.in.Degrees.C'")
      message(paste("Number of missing values in Dew.point.temperature.in.Degrees.C: ", sum(is.na(climatedata$Solar.radiation.in.MJ.per.sqm.per.day))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Solar.radiation.in.MJ.per.sqm.per.day))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Solar.radiation.in.MJ.per.sqm.per.day)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Solar.radiation.in.MJ.per.sqm.per.day exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    Rs.temp <- zoo(climatedata$Solar.radiation.in.MJ.per.sqm.per.day, Date.subdaily)
    Rs <- aggregate(Rs.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message("Interpolation is performed to replace data with error, users should be aware of the associated risk")
    for (m in 0:11) {
      Rs[as.POSIXlt(time(Rs))$mon==m & is.na(Rs)] = mean(Rs[as.POSIXlt(time(Rs))$mon==m & !is.na(Rs)])
    }
  } else {
      Rs <- NULL
  }
    # Check if data 'n' exists
    if ("n" %in% (colnames(climatedata))) {
      if ("TRUE" %in% (is.na(climatedata$n))) {
        message("Warning: missing values in 'n' (daily sunshine hours)")
        message(paste("Number of missing values in n: ", sum(is.na(climatedata$n))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$n))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$n)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of n exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      n.temp <- zoo(climatedata$n, Date.subdaily)  
      for (m in 0:11) {
        n.temp[as.POSIXlt(time(n.temp))$mon==m &  as.numeric(is.na(n.temp))] = mean(n.temp[as.POSIXlt(time(n.temp))$mon==m &  as.numeric(!is.na(n.temp))]) # Changing to monthly mean 
      }
      n <- aggregate(n.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    } else {
      n <- NULL
    } 
  
    # Check if data 'cloud.cover.in.oktas' exists
    if ("cloud.cover.in.oktas" %in% (colnames(climatedata))) {
      if ("TRUE" %in% (is.na(climatedata$cloud.cover.in.oktas))) {
        message("Warning: missing values in 'cloud.cover.in.oktas' (daily cloud cover)")
        message(paste("Number of missing values in cloud.cover.in.oktas: ", sum(is.na(climatedata$cloud.cover.in.oktas))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$cloud.cover.in.oktas))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$cloud.cover.in.oktas)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of cloud.cover.in.oktas exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      C0.temp <- zoo(climatedata$cloud.cover.in.oktas, Date.subdaily) 
      for (m in 0:11) {
        C0.temp[as.POSIXlt(time(C0.temp))$mon==m &  as.numeric(is.na(C0.temp))] = mean(C0.temp[as.POSIXlt(time(C0.temp))$mon==m &  as.numeric(!is.na(C0.temp))]) # Changing to monthly mean 
      }
      C0 <- aggregate(C0.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
      n <- constants$a_0 + constants$b_0*C0 + constants$c_0*C0^2 + constants$d_0*C0^3 # calculation of sunshine hours (h) based on cloud cover (oktas) (S3.10)
    } 
  
  
    # Check if data 'precipitation.in.mm' exists
    if ("precipitation.in.mm" %in% (colnames(climatedata))) {
      if ("TRUE" %in% (is.na(climatedata$precipitation.in.mm))) {
        message("Warning: missing values in 'precipitation.in.mm' (daily precipitation)")
        message(paste("Number of missing values in precipitation.in.mm: ", sum(is.na(climatedata$precipitation.in.mm))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$precipitation.in.mm))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$precipitation.in.mm)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of precipitation.in.mm exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      P.temp <- zoo(climatedata$precipitation.in.mm, Date.subdaily) 
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
      Cd.temp <- zoo(Cd.temp, Date.daily)
      Cd <- Precip
      for (m in 1:length(Cd)) {
        Cd[as.yearmon(time(Cd)) == as.yearmon(time(Cd.temp))[m]] <- Cd.temp[m]
      }
    } else {
      Cd <- NULL
    } 

 
  # Check if data 'precipitation.in.mm' exist
  if ("precipitation.in.mm" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$precipitation.in.mm))) {
      message("Warning: missing values in 'precipitation.in.mm' (daily precipitation)")
      message(paste("Number of missing values in precipitation.in.mm: ", sum(is.na(climatedata$precipitation.in.mm))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$precipitation.in.mm))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$precipitation.in.mm)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of precipitation.in.mm exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    P.temp <- zoo(climatedata$precipitation.in.mm, Date.subdaily) 
    for (m in 0:11) {
      P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(is.na(P.temp))] = mean(P.temp[as.POSIXlt(time(P.temp))$mon==m &  as.numeric(!is.na(P.temp))]) # Changing to monthly mean 
    }
    Precip <- aggregate(P.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
  } else {
    Precip <- NULL
  }

  # Check if data 'Class.A.pan.evaporation.in.mm' exist
  if ("Class.A.pan.evaporation.in.mm" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Class.A.pan.evaporation.in.mm))) {
        message("Warning: missing values in 'Class.A.pan.evaporation.in.mm'")
        message(paste("Number of missing values in Class.A.pan.evaporation.in.mm: ", sum(is.na(climatedata$Class.A.pan.evaporation.in.mm))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Class.A.pan.evaporation.in.mm))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$Class.A.pan.evaporation.in.mm)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of Class.A.pan.evaporation.in.mm exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      Epan.temp <- zoo(climatedata$Class.A.pan.evaporation.in.mm, Date.subdaily) 
      for (m in 0:11) {
        Epan.temp[as.POSIXlt(time(Epan.temp))$mon==m &  as.numeric(is.na(Epan.temp))] = mean(Epan.temp[as.POSIXlt(time(Epan.temp))$mon==m &  as.numeric(!is.na(Epan.temp))]) # Changing to monthly mean 
      }
      Epan <- aggregate(P.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean) 
  } else {
    Epan <- NULL
  }
  
  # Check if data 'RHmax' and 'RHmin' exist
  if ("RHmax" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$RHmax))) {
      message("Warning: missing values in 'RHmax' (daily maximum temperature)")
      message(paste("Number of missing values in RHmax: ", sum(is.na(climatedata$RHmax))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmax))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$RHmax)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of RHmax exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    RHmax.temp <- zoo(climatedata$RHmax, Date.subdaily)
    RHmax <- aggregate(RHmax.temp, as.Date(Date.subdaily, "%d/%m/%y") ,FUN = mean) 
    for (m in 0:11) {
      RHmax[as.POSIXlt(time(RHmax))$mon==m & is.na(RHmax)] = mean(RHmax[as.POSIXlt(time(RHmax))$mon==m & !is.na(RHmax)])
    }
    message(paste("Number of days increments when RHmax has errors: ", sum(RHmax>100)))
    message("'RHmax' with values > 100% has been adjusted to 100%, users should be aware of the associated risk")
    RHmax[RHmax>100] = 100
  } else {
    if ("Relative.humidity.in.percentage.." %in% (colnames(climatedata))) {
      message("Warning: missing data of 'RHmax'(daily maximum relative humidity), calculated from subdaily 'Relative.humidity.in.percentage..'")
      if ("TRUE" %in% (is.na(climatedata$Relative.humidity.in.percentage..))) {
        message("Warning: missing values in 'Relative.humidity.in.percentage..'")
        message(paste("Number of missing values in Relative.humidity.in.percentage..: ", sum(is.na(climatedata$Relative.humidity.in.percentage..))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Relative.humidity.in.percentage..))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$Relative.humidity.in.percentage..)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of Relative.humidity.in.percentage.. exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
    RH.temp <- zoo(climatedata$Relative.humidity.in.percentage.., Date.subdaily)
    RHmax <- aggregate(RH.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = max)
    for (m in 0:11) {
      RHmax[as.POSIXlt(time(RHmax))$mon==m & is.na(RHmax)] = mean(RHmax[as.POSIXlt(time(RHmax))$mon==m & !is.na(RHmax)])
    }
    message(paste("Number of days increments when RHmax has errors: ", sum(RHmax>100)))
    message("'RHmax' with values > 100% has been adjusted to 100%, users should be aware of the associated risk")
    RHmax[RHmax>100] = 100
      } else {
      RHmax <- NULL
    }
  }
  
  if ("RHmin" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$RHmin))) {
      message("Warning: missing values in 'RHmin' (daily minimum temperature)")
      message(paste("Number of missing values in RHmin: ", sum(is.na(climatedata$RHmin))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$RHmin))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$RHmin)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of RHmin exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    RHmin.temp <- zoo(climatedata$RHmin, Date.subdaily)
    RHmin <- aggregate(RHmin.temp, as.Date(Date.subdaily, "%d/%m/%y") ,FUN = mean) 
    for (m in 0:11) {
      RHmin[as.POSIXlt(time(RHmin))$mon==m & is.na(RHmin)] = mean(RHmin[as.POSIXlt(time(RHmin))$mon==m & !is.na(RHmin)])
    }
    message(paste("Number of days increments when RHmin has errors: ", sum(RHmin>RHmax)))
    message("'RHmax' with values > 'RHmax' has been adjusted to '0.9999*RHmax', users should be aware of the associated risk")
    RHmin[RHmin>=RHmax] = 0.9999*RHmax[RHmin>=RHmax] # setting upper bound on minimum RH to ensure that ea < es for all es,ea. 
    } else {
    if ("Relative.humidity.in.percentage.." %in% (colnames(climatedata))) {
      message("Warning: missing data of 'RHmin'(daily minimum relative humidity), calculated from subdaily 'Relative.humidity.in.percentage..'")
      if ("TRUE" %in% (is.na(climatedata$Relative.humidity.in.percentage..))) {
        message("Warning: missing values in 'Relative.humidity.in.percentage..'")
        message(paste("Number of missing values in Relative.humidity.in.percentage..: ", sum(is.na(climatedata$Relative.humidity.in.percentage..))))
        message(paste("% missing data: ", signif(sum(is.na(climatedata$Relative.humidity.in.percentage..))/nrow(climatedata) * 100, digits = -3), "%"))
        if (sum(is.na(climatedata$Relative.humidity.in.percentage..)) >= stopmissing/100 * nrow(climatedata)) {
          stop("missing data of Relative.humidity.in.percentage.. exceeds ", stopmissing, "%, please use high quality data for calculation")
        } else {
          message("Interpolation is performed for missing data, users should be aware of the associated risk")
        }
      }
      RH.temp <- zoo(climatedata$Relative.humidity.in.percentage.., Date.subdaily)
      RHmin <- aggregate(RH.temp, as.Date(Date.subdaily, "%d/%m/%y"), FUN = min)
      for (m in 0:11) {
        RHmin[as.POSIXlt(time(RHmin))$mon==m & is.na(RHmin)] = mean(RHmin[as.POSIXlt(time(RHmin))$mon==m & !is.na(RHmin)])
      }
      message(paste("Number of days increments when RHmin has errors: ", sum(RHmin>RHmax)))
      message("'RHmax' with values > 'RHmax' has been adjusted to '0.9999*RHmax', users should be aware of the associated risk")
      RHmin[RHmin>=RHmax] = 0.9999*RHmax[RHmin>=RHmax] # setting upper bound on minimum RH to ensure that ea < es for all es,ea. 
      } else {
      RHmin <- NULL
    }
  }
  
  # Check if data 'Dew.point.temperature.in.Degrees.C' exists
  if ("Dew.point.temperature.in.Degrees.C" %in% (colnames(climatedata))) {  
    if ("TRUE" %in% (is.na(climatedata$Dew.point.temperature.in.Degrees.C))) {
      message("Warning: missing values in 'Dew.point.temperature.in.Degrees.C'")
      message(paste("Number of missing values in Dew.point.temperature.in.Degrees.C: ", sum(is.na(climatedata$Dew.point.temperature.in.Degrees.C))))
      message(paste("% missing data: ", signif(sum(is.na(climatedata$Dew.point.temperature.in.Degrees.C))/nrow(climatedata) * 100, digits = -3), "%"))
      if (sum(is.na(climatedata$Dew.point.temperature.in.Degrees.C)) >= stopmissing/100 * nrow(climatedata)) {
        stop("missing data of Dew.point.temperature.in.Degrees.C exceeds ", stopmissing, "%, please use high quality data for calculation")
      } else {
        message("Interpolation is performed for missing data, users should be aware of the associated risk")
      }
    }
    Tdew.temp <- zoo(climatedata$Dew.point.temperature.in.Degrees.C, Date.subdaily)
    Tdew <- aggregate(Tdew.temp, as.Date(Date.subdaily, "%d/%m/%y"),mean)
    message(paste("Number of days increments when Tdew has errors: ", sum(Tdew>100)))
    message("Interpolation is performed to replace data with error, users should be aware of the associated risk")
    for (m in 0:11) {
      Tdew[as.POSIXlt(time(Tdew))$mon==m & is.na(Tdew)] = mean(Tdew[as.POSIXlt(time(Tdew))$mon==m & !is.na(RHmin)])
      Tdew[as.POSIXlt(time(Tdew))$mon==m & as.numeric(Tdew)>100] = mean(Tdew[as.POSIXlt(time(Tdew))$mon==m & as.numeric(Tdew)<100])
    }
    } else {
      if ("Vapour.pressure.in.hPa" %in% (colnames(climatedata))) {
        message("Warning: missing data of 'Dew.point.temperature.in.Degrees.C', calculated from 'Vapour.pressure.in.hPa'")
        Tdew <- NULL
        if ("TRUE" %in% (is.na(climatedata$Vapour.pressure.in.hPa))) {
          message("Warning: missing values in 'Vapour.pressure.in.hPa'")
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

PlotEvapotranspiration <- function(results, OBS, OBSplot) { # plot estimations and observations
  # Plots - Aggregation
  par(ask=FALSE)
  plot.new()
  if (!is.null(results$PET.Daily)) {
    plot(results$PET.Daily, main = paste("Daily", results$PET_formulation, results$PET_type),  ylim = c(0,20), xlab = "Year", ylab = list(c(results$PET_type, "mm/day")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        if (!is.null(OBS$E_obs.Daily)) {
          lines(OBS$E_obs.Daily, main = "Observed evaporation", type = "o", pch = ".", ylim = c(0,400), col = "RED", xlab = "Year", ylab = "Observed evaporation mm/day")
          legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
        }
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
  } else {
    warning("No daily results obtained from ", results$PET_formulation)
  }
  par(ask=TRUE)
  
  plot(results$PET.Monthly, main = paste("Monthly", results$PET_formulation, results$PET_type), type = "o", pch = ".", ylim = c(0,400), xlab = "Year", ylab = list(c(results$PET_type, "mm/month")))
  if (OBSplot == TRUE) {
    if (!is.null(OBS)) {
      lines(OBS$E_obs.Monthly, main = "Observed evaporation", type = "o", pch = ".", ylim = c(0,400), col = "RED", xlab = "Year", ylab = "Observed evaporation mm/month")
      legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
    } else {
      warning("No observed data available for plotting, only the estimated values are plotted")
    }
  }
  plot(results$PET.Annual, main = paste("Annual", results$PET_formulation, results$PET_type), type = "o", xlab = "Year", ylim = c(0,3000), ylab = list(c(results$PET_type, "mm/year")))
  if (OBSplot == TRUE) {
    if (!is.null(OBS)) {
      lines(OBS$E_obs.Annual, main = "Observed evaporation", type = "o", col = "RED", ylim = c(0,3000), xlab = "Year", ylab = "Observed evaporation mm/year")
      legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
    } else {
      warning("No observed data available for plotting, only the estimated values are plotted")
    }
  }
  
  # Plots - Average
  plot(results$PET.MonthlyAve~unique(1:12), main = paste("Monthly average", results$PET_formulation, results$PET_type), type = "o", ylim = c(0,10), xlab = "Month", ylab = list(c("Monthly average", results$PET_type, "mm/day")))
  if (OBSplot == TRUE) {
    if (!is.null(OBS)) {
      lines(OBS$E_obs.MonthlyAve~unique(1:12), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Month", ylab = "Monthly average observed evaporation mm")
      legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
    } else {
      warning("No observed data available for plotting, only the estimated values are plotted")
    }
  }
  plot(results$PET.AnnualAve~unique(as.POSIXlt(data$Date.daily)$year + 1900), main = paste("Annual average", results$PET_formulation, results$PET_type), type = "o", ylim = c(0,10), xlab = "Year", ylab = list(c("Annual average", results$PET_type, "mm/day")))
  if (OBSplot == TRUE) {
    if (!is.null(OBS)) {
      lines(OBS$E_obs.AnnualAve~unique(as.POSIXlt(OBS$Date.OBS)$year + 1900), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Year", ylab = "Annual average observed evaporation mm")
      legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
    } else {
      warning("No observed data available for plotting, only the estimated values are plotted")
    }
  }
  par(ask=FALSE)
  paste("Completed plotting for results calculated by", results$PET_formulation, "formulation")
}

  #-------------------------------------------------------------------------------------

EvapotranspirationForcings <- function(data, results, forcing) {
  plot.new()
  # define forcing names
  Fnames <- list(v=c("Tmax","Tmin","u2","uz","Rs","n","Precip","Epan","RHmax","RHmin","Tdew"), 
              n=c("maximum temperature","minimum temperature", "average wind speed at 2m", "average wind speed", "daily solar radiation", "daily sunshine hours", 
                 "precipitation", "Class-A pan evaporation", "maximum relative humidity", "minimum relative humidity", "average dew point temeprature"),
              u=c("degree Celcius","degree Celcius", "m/s", "m/s", "MJ.per.sqm.per.day", "hours", "mm/day", "mm/day", "%", "%", "degree Celcius"))
  
  # check forcing value
  if (forcing %in% Fnames$v) {
    par(mfrow=c(2,2))
    
    # Aggregated PET vs forcing
    
    # Plot daily relationship
    if (!is.null(results$PET.Daily)) {
      interval <- "daily"
      if (!is.null(data[[forcing]])) {
        plot(results$PET.Daily~data[[forcing]], main = results$PET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$PET_type, "mm/day")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no daily ET estimation to plot")
    }
  
    # Plot monthly relationship
    if (!is.null(results$PET.Monthly)) {
      interval <- "monthly"
      if (!is.null(data[[forcing]])) {
        if (forcing != "Precip") {
          Fmonthly <- aggregate(data[[forcing]], as.yearmon(data$Date.daily, "%m/%y"), mean)
        } else {
          Fmonthly <- aggregate(data[[forcing]], as.yearmon(data$Date.daily, "%m/%y"), sum)
        }
        plot(results$PET.Monthly~Fmonthly, main = results$PET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$PET_type, "mm/month")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no monthly ET estimation to plot")
    }
    
    
    # Plot annual relationship
    if (!is.null(results$PET.Annual)) {
      interval <- "annual"
      if (!is.null(data[[forcing]])) {
        if (forcing != "Precip") {
          Fannual <- aggregate(Fmonthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), mean)
        } else {
          Fannual <- aggregate(Fmonthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), sum)
        }
        plot(results$PET.Annual~Fannual, main = results$PET_formulation, xlab = paste(interval, Fnames$n[Fnames$v == forcing], Fnames$u[Fnames$v == forcing]), ylab = list(c(results$PET_type, "mm/year")))
      } else {
        plot(1, type="n", axes=F, xlab="", ylab="")
        text(1, paste("no", forcing, "data to plot"))
      }
    } else {
      plot(1, type="n", axes=F, xlab="", ylab="")
      text(1, "no annual ET estimation to plot")
    }
  
    
  } else {
    stop("Please select forcing from 'Tmax','Tmin','u2','uz','Rs','n','Precip','Epan','RHmax','RHmin','Tdew'") # forcing is not defined
  }
  
}

  #-------------------------------------------------------------------------------------

Evapotranspiration <- function(data, ...) UseMethod("Evapotranspiration")

  #-------------------------------------------------------------------------------------
  
Evapotranspiration.Penman <- function(data, constants, solar, wind, windfunction_ver, alpha = 0.08, z0 = 0.001, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (daily average temperature)")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
    }
    if (is.null(data$uz) & is.null(data$u2)) {
      stop("fail to obtain data of wind speed")
    }
  }

  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  } 
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  
  # check user-input albedo
  if (wind == "yes") {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    }
    if (is.na(as.numeric(z0))) {
      stop("Please use a numeric value for the z0 (roughness height)")
    }  
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Penman
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
    if (solar == "data") {
      R_s <- data$Rs
    } else if (solar!="monthly precipitation") {
      # calculate R_s from sunshine hours - data or estimation using cloudness
      R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
    } else {
      # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
      R_s <- (0.85 - 0.047*data$Cd)*R_a 
    }
    
    if (wind == "yes") {
      # Wind speed at 2 meters
      if (is.null(data$u2)) {
        u2 <- data$uz * log(2/z0) / log(constants$z/z0) # Equation S4.4
      } else {
        u2 <- data$u2
      }
      
      # Saturated vapour pressure
      vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
      vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
      vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
      
      # Vapour pressure
      vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
      
      R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
      R_ns <- (1 - alpha) * R_s # net incoming shortwave radiation - water or other evaporative surface with specified Albedo (S3.2)
      
      R_n = R_ns - R_nl # net radiation (S3.1)
      if (windfunction_ver == "1948") {
        f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
      } else if (windfunction_ver == "1956") {
        f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
      } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
        stop("Please select the version of wind function (1948 or 1956)")
      }
      
      Ea = f_u * (vas - vabar) # (S4.2)
      
      Epenman.Daily <-  delta / (delta +  gamma) * (R_n / constants$lambda) + gamma  / (delta + gamma) * Ea # Penman open-water evaporation (S4.1)
    } else {
      # mean relative humidity
      RHmean <- (data$RHmax + data$RHmin) / 2 
      
      Epenman.Daily <-  0.047 * R_s * sqrt(Ta + 9.5) - 2.4 * (R_s/R_a)^2 + 0.09 * (Ta + 20) * (1 - RHmean/100) # Penman open-water evaporation without wind data by Valiantzas (2006) (S4.12)
    }
   
   PET.Daily <- Epenman.Daily
   PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
   PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
   PET.MonthlyAve <- PET.AnnualAve <- NULL
   for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
     i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
     PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
   }
   for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
     i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
     PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
   }

   # Generate summary message for results  
   PET_formulation <- "Penman" 
   if (wind == "no") {
     PET_type <- "Open-water Evaporation"
     Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
   } else {
     if (alpha != 0.08) {
       PET_type <- "Potential Evaporation"
       Surface <- paste("user-defined, albedo =", alpha, "; roughness height =", z0, "m")
     } else if (alpha == 0.08) {
       PET_type <- "Open-water Evaporation"
       Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
     }
   }
   
   if (solar == "data") {
     message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"    
   } else if (solar == "sunshine hours") {
     message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
   } else if (solar == "cloud") {
     message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
   } else {
     message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
   }
   
   if (wind == "yes") {
     if (windfunction_ver == "1948") {
       message2 <- "Wind data has been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
     } else if (windfunction_ver == "1956") {
       message2 <- "Wind data has been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
     } 
   } else {
     message2 <- "Alternative calculation for Penman evaporation without wind data has been performed"
   }
  
   message(PET_formulation, " ", PET_type)
   message("Evaporative surface: ", Surface)
   message(message1)
   message(message2)
  
   results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message2=message2)
   class(results) <- funname
   return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.PenmanMonteith <- function(data, constants, solar, wind, crop, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
  stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("fail to obtain data of wind speed")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  } 
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  # check user-input crop type and specify albedo
  if (wind == "yes") {
    if (crop != "short" & crop != "tall") {
      stop("Please enter 'short' or 'tall' for the desired reference crop type")
    } else {
      alpha <- 0.23 # albedo for both short and tall crop
      if (crop == "short") {
        z0 <- 0.02 # roughness height for short grass
      } else {
        z0 <- 0.1 # roughness height for tall grass
      }
    }
  } else {
    z0 <- 0.02 # roughness height for short grass
    alpha <- 0.25 # semi-desert short grass - will not be used for calculation - just informative
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Penman-Monteith Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  if (wind == "yes") {
    # Wind speed
    if (is.null(data$u2)) {
      u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
    } else {
      u2 <- data$u2
    }
    
    if (crop == "short") {
      r_s <- 70 # will not be used for calculation - just informative
      CH <- 0.12 # will not be used for calculation - just informative
      ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) + gamma * 900 * u2 * (vas - vabar)/(Ta + 273)) / (delta + gamma * (1 + 0.34*u2)) # FAO-56 reference crop evapotranspiration from short grass (S5.18)
    } else {
      r_s <- 45 # will not be used for calculation - just informative
      CH <- 0.50 # will not be used for calculation - just informative
      ET_RC.Daily <- (0.408 * delta * (R_ng - constants$G) + gamma * 1600 * u2 * (vas - vabar)/(Ta + 273)) / (delta + gamma * (1 + 0.38*u2)) # ASCE-EWRI standardised Penman-Monteith for long grass (S5.19)
    }
  } else {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    R_s.Monthly <- aggregate(R_s, as.yearmon(data$Date.daily, "%m/%y"),mean)
    R_a.Monthly <- aggregate(R_a, as.yearmon(data$Date.daily, "%m/%y"),mean)
    Ta.Monthly <- aggregate(Ta, as.yearmon(data$Date.daily, "%m/%y"),mean)
    RHmean.Monthly <- aggregate(RHmean, as.yearmon(data$Date.daily, "%m/%y"),mean)
    ET_RC.Daily <- 0.038 * R_s.Monthly * sqrt(Ta.Monthly + 9.5) - 2.4 * (R_s.Monthly/R_a.Monthly)^2 + 0.075 * (Ta.Monthly + 20) * (1 - RHmean.Monthly/100) # Reference crop evapotranspiration without wind data by Valiantzas (2006) (S5.21)
  }
  
  PET.Daily <- ET_RC.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (wind == "no") {
    PET_formulation <- "Penman-Monteith (without wind data)"
    PET_type <- "Reference Crop Evapotranspiration"
    Surface <- paste("short grass, albedo =", alpha, "; roughness height =", z0, "m")
  } else {
    if (crop == "short") {
      PET_formulation <- "Penman-Monteith FAO56"
      PET_type <- "Reference Crop Evapotranspiration"
      Surface <- paste("FAO-56 hypothetical short grass, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    } else {
      PET_formulation <- "Penman-Monteith ASCE-EWRI Standardised"
      PET_type <- "Reference Crop Evapotranspiration"
      Surface <- paste("ASCE-EWRI hypothetical tall grass, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    }
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    message2 <- "Wind data has been used for calculating the reference crop evapotranspiration"
  } else {
    message2 <- "Alternative calculation for reference crop evapotranspiration without wind data has been performed"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message2=message2)
  class(results) <- funname
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.MattShuttleworth <- function(data, constants, solar, alpha, r_s, CH, ...) { 
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature), 'vas' (daily saturated vapour pressure) and 'vabar' (mean daily vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("fail to obtain data of wind speed")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # check user-input albedo, surface resistance and crop height
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (is.na(as.numeric(r_s))) {
    stop("Please use a numeric value for the r_s (surface resisitance) in sm^-1")
  }
  if (is.na(as.numeric(CH))) {
    stop("Please use a numeric value for the CH (crop height) in m")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Matt-Shuttleworth Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  r_clim <- 86400 * constants$Roua * constants$Ca * (vas - vabar) / delta * R_ng # clinmatological resistance (s*m^-1) (S5.34)
  r_clim[r_clim == 0] <- 0.1 # correction for r_clim = 0
  u2[u2 == 0] <- 0.1 # correction for u2 = 0
  VPD50toVPD2 <- (302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) + 1/r_clim * ((302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) * (208 / u2) - (302 / u2)) # ratio of vapour pressure deficits at 50m to vapour pressure deficits at 2m heights (S5.35)
  r_c50 <- 1 / ((0.41)^2) * log((50 - 0.67 * CH) / (0.123 * CH)) * log((50 - 0.67 * CH) / (0.0123 * CH)) * log((2 - 0.08) / 0.0148) / log((50 - 0.08) / 0.0148) # aerodynamic coefficient (s*m^-1) (S5.36)
  
  E_Tc.Daily <- 1/constants$lambda * (delta * R_ng + (constants$Roua * constants$Ca * u2 * (vas - vabar)) / r_c50) * VPD50toVPD2 / (delta + gamma * (1 + r_s * u2 / r_c50)) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)
  
  PET.Daily <- E_Tc.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Matt-Shuttleworth"
  PET_type <- "Reference Crop Evapotranspiration"
  Surface <- paste("user-defined, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.PriestleyTaylor <- function(data, constants, solar, alpha, ...) {  
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature) and 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 

  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7

  # Calculations from data and constants for Matt-Shuttleworth Reference Crop
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  
  E_PT.Daily <- constants$alphaPT * (delta/(delta + gamma) * R_ng / constants$lambda - constants$G / constants$lambda) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)
  
  PET.Daily <- E_PT.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Priestley-Taylor"
  PET_type <- "Potential Evaporation"
  if (alpha != 0.08) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else if (alpha == 0.08) {
    Surface <- paste("water, albedo =", alpha)
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.Penpan <- function(data, constants, solar, alpha, overest, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature), 'vas' (daily saturated vapour pressure) and 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("fail to obtain data of wind speed")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Penpan
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) * (constants$lat)^2 # pan radiation factor (S6.6)
  f_dir <- -0.11 + 1.31 * R_s / R_a # fraction of R_S that os direct (S6.5)
  R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) * R_s # total shortwave radiation received (S6.4)
  R_npan <-(1 - constants$alphaA) * R_span - R_nl # net radiation at the pan (S6.3)
  f_pan_u <-1.201 + 1.621 * u2 # (S6.2)
  
  Epenpan.Daily <- delta / (delta + constants$ap * gamma) * R_npan / constants$lambda + constants$ap * gamma / (delta + constants$ap * gamma) * f_pan_u * (vas - vabar) # Penpan estimation of Class-A pan evaporation (S6.1)
  if (overest == TRUE) {
    Epenpan.Daily <- Epenpan.Daily / 1.078 # adjusted by 1.078 to balance overestimation observed by McMahon 
  }
  
    
  PET.Daily <- Epenpan.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Penpan"
  PET_type <- "Class-A Pan Evaporation"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.BrutsaertStrickler <- function(data, constants, solar, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature), 'vabar' (mean daily actual vapour pressure) and 'vas' (daily saturated vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("fail to obtain data of wind speed")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # update alphaPT according to Brutsaert and Strickler (1979)
  constants$alphaPT <- 1.28
  
  # Calculations from data and constants for Brutsaert and Strickler actual evapotranspiration
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
  f_u2 <- 2.626 + 1.381 * u2 # Penman's wind function adopted with Priestley and Tayloy constant alphaPT by Brutsaert and Strickler (1979) (S8.3)
  
  ET_BS_Act.Daily <- (2 * constants$alphaPT - 1) * (delta / (delta + gamma)) * R_ng / constants$lambda - gamma / (delta + gamma) * f_u2 * (vas - vabar) # Brutsaert and Strickler actual areal evapotranspiration (mm.day^-1) Brutsaert and Strickler (1979) (S8.2)
  
  PET.Daily <- ET_BS_Act.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Brutsaert-Strickler"
  PET_type <- "Actual Areal Evapotranspiration"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.GrangerGray <- function(data, constants, solar, windfunction_ver, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature), 'vabar' (mean daily actual vapour pressure) and 'vas' (daily saturated vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("fail to obtain data of wind speed")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # check user-input albedo
  if (is.na(as.numeric(alpha))) {
    stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
  }
  if (!is.na(as.numeric(alpha))) {
    if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
      stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
    }
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Granger and Gray actual evapotranspiration
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
  # For short grass
  R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
  R_ng <- R_nsg - R_nl # net radiation (S3.1)
 
  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  if (windfunction_ver == "1948") {
    f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
  } else if (windfunction_ver == "1956") {
    f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
  } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
    stop("Please select the version of wind function (1948 or 1956)")
  }
  Ea = f_u * (vas - vabar) # (S4.2)
  D_p <- Ea / (Ea + (R_ng - constants$G) / constants$lambda) # dimensionless relative drying power (S8.6)
  G_g <- 1 / (0.793 + 0.20 * exp(4.902 * D_p)) + 0.006 * D_p # dimensionless evaporation parameter (S8.5)
  ET_GG_Act.Daily <- delta * G_g / (delta * G_g + gamma) * (R_ng - constants$G) / constants$lambda + gamma * G_g / (delta * G_g + gamma) * Ea # Granger and Gray actual areal evapotranspiration (mm.day^-1) Granger and Gray (1989) (S8.4)
  
  PET.Daily <- ET_GG_Act.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Granger-Gray"
  PET_type <- "Actual Areal Evapotranspiration"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (windfunction_ver == "1948") {
    message2 <- "Wind data has been used for the calculation of the drying power of air, using Penman 1948 wind function."
  } else if (windfunction_ver == "1956") {
    message2 <- "Wind data has been used for the calculation of the drying power of air, using Penman 1956 wind function."
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)

  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message2=message2)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.SzilagyiJozsa <- function(data, constants, solar, wind, windfunction_ver, alpha, z0, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature) and 'vabar' (mean daily actual vapour pressure)")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("fail to obtain data of wind speed")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  if (wind != "yes" & wind != "no") {
    stop("Please choose if actual data will be used for wind speed from wind = 'yes' and wind = 'no'")
  }
  
  # check user-input albedo
  if (wind == "yes") {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    }
    if (is.na(as.numeric(z0))) {
      stop("Please use a numeric value for the z0 (roughness height)")
    }  
  }

  # update alphaPT according to Szilagyi and Jozsa (2008)
  constants$alphaPT <- 1.31
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Saturated vapour pressure
  vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
  vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
  vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
  
  # Vapour pressure
  vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
  
  # Calculations from data and constants for Penman
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  if (wind == "yes") {
    # Wind speed at 2 meters
    if (is.null(data$u2)) {
      u2 <- data$uz * log(2/z0) / log(constants$z/z0) # Equation S4.4
    } else {
      u2 <- data$u2
    }
    
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For vegetated surface
    R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
    R_ng = R_nsg - R_nl # net radiation (S3.1)
    if (windfunction_ver == "1948") {
      f_u = 2.626 + 1.381 * u2 # wind function Penman 1948 (S4.11)
    } else if (windfunction_ver == "1956") {
      f_u = 1.313 + 1.381 * u2 # wind function Penman 1956 (S4.3)
    } else if (windfunction_ver != "1948" & windfunction_ver != "1956") {
      stop("Please select the version of wind function (1948 or 1956)")
    }
    Ea = f_u * (vas - vabar) # (S4.2)
    
    Epenman.Daily <-  delta / (delta +  gamma) * (R_ng / constants$lambda) + gamma  / (delta + gamma) * Ea # Penman open-water evaporation (S4.1)
  } else {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    Epenman.Daily <-  0.047 * R_s * sqrt(Ta + 9.5) - 2.4 * (R_s/R_a)^2 + 0.09 * (Ta + 20) * (1 - RHmean/100) # Penman open-water evaporation without wind data by Valiantzas (2006) (S4.12)
  }
  
  # Iteration for equilibrium temperature T_e
  T_e <- Ta
  for (i in 1:99999) {
    v_e <- 0.6108 * exp(17.27 * T_e/(T_e + 237.3)) # saturated vapour pressure at T_e (S2.5)
    T_enew <- Ta - 1 / gamma * (1 - R_ng / (constants$lambda * Epenman.Daily)) * (v_e - vabar) # rearranged from S8.8
    deltaT_e <- na.omit(T_enew - T_e)
    maxdeltaT_e <- abs(max(deltaT_e))
    T_e <- T_enew
    if (maxdeltaT_e < 0.01) break
  }
  deltaT_e <- 4098 * (0.6108 * exp((17.27 * T_e)/(T_e+237.3))) / ((T_e + 237.3)^2)  # slope of vapour pressure curve (S2.4)
  E_PT_T_e <- constants$alphaPT * (deltaT_e / (deltaT_e + gamma) * R_ng / constants$lambda) # Priestley-Taylor evapotranspiration at T_e
  E_SJ_Act.Daily <- 2 * E_PT_T_e - Epenman.Daily # actual evapotranspiration by Szilagyi and Jozsa (2008) (S8.7)
  
  PET.Daily <- E_SJ_Act.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Szilagyi-Jozsa"
  PET_type <- "Actual Evapotranspiration"
  Surface <- paste("user-defined, albedo =", alpha, "; roughness height", z0, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    if (windfunction_ver == "1948") {
      message2 <- "Wind data has been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
    } else if (windfunction_ver == "1956") {
      message2 <- "Wind data has been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
    } 
  } else {
    message2 <- "Alternative calculation for Penman evaporation without wind data has been performed"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message2=message2)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.Makkink <- function(data, constants, solar, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Calculations from data and constants for Makkink
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  Emakkink.Daily <- 0.61 * (delta / (delta + gamma) * R_s/2.45) - 0.12  # potential evapotranspiration by Bruin (1981) (S9.6)
  
  PET.Daily <- Emakkink.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Makkink"
  PET_type <- "Potential Evaporation"

  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: open-water")
  message(message1)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.BalneyCriddle <- function(data, constants, solar, height, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  if (solar == "sunshine hours" & is.null(data$n)) { # sunshine hour data is required
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud") {
    if (is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
      stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("fail to obtain data of wind speed")
    }
    if (is.null(data$RHmin)) { 
      stop("fail to obtain 'RHmin' (minimum daily relative humidity)")
    } 
  } 
  if (solar == "data" | solar == "monthly precipitation") {
    stop("Only 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
  }
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Calculations from data and constants for Balney and Criddle
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values

  # Wind speed
  if (is.null(data$u2)) {
    u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
  } else {
    u2 <- data$u2
  }
  
  bvar <- constants$e0 + constants$e1 * data$RHmin + constants$e2 * data$n/N + constants$e3 * u2 + constants$e4 * data$RHmin * data$n/N + constants$e5 * data$RHmin * u2 # undefined working variable (Allena and Pruitt, 1986; Shuttleworth, 1992) (S9.8)
  N.annual <- ave(N, format(time(N), "%y"), FUN = sum) # Annual sum of maximum sunshine hours
  p_y <- 100 * data$n/N.annual # percentage of actual daytime hours for the day comparing to the annual sum of maximum sunshine hours

  
  ET_BC.Daily <- (0.0043 * data$RHmin - data$n/N - 1.41) + bvar * p_y * (0.46 * Ta +8.13) # Balney-Criddle Reference Crop evapotranspiration (mm.day^-1) (S9.7)
  
  PET.Daily <- ET_BC.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  if (height == TRUE) {
    ET_BC.Daily = ET_BC.Daily * (1 + 0.1 * constants$Elev/1000) # with adjustment for site elevation by Allen and Pruitt (1986) (S9.9) 
  }
  
  # Generate summary message for results
  PET_formulation <- "Balney-Criddle"
  PET_type <- "Reference Crop Evapotranspiration"

  if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (height == TRUE) {
    message3 <- "Height adjustment has been applied to calculated Balney-Criddle reference crop evapotranspiration"
  } else {
    message3 <- "No height adjustment has been applied to calculated Balney-Criddle reference crop evapotranspiration"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message3)
    
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message3=message3)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.Truc <- function(data, constants, solar, humid, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("fail to obtain data of 'Rs' (daily solar radiation)")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("fail to obtain data of 'n' (daily sunshine hours)")
  } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
    stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
  } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
  }

  if (humid == TRUE & (is.null(data$RHmax)|is.null(data$RHmin))) { # for adjustment for non-humid conditions
    stop("fail to obtain 'RHmean' (average daily relative humidity)")
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Truc
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar!="monthly precipitation") {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  } else {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  }
  
  ET_truc.Daily <- 0.013 * (23.88 * R_s + 50) * Ta / (Ta + 15) # reference crop evapotranspiration by Truc (1961) (S9.10)
  
  if (humid == TRUE) {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    ET_truc.Daily[RHmean < 50] <- 0.013 * (23.88 * R_s + 50) * Ta[RHmean < 50] / (Ta[RHmean < 50] + 15) * (1 + (50 - RHmean[RHmean < 50]) / 70) # Truc reference crop evapotranspiration adjusted for non-humid conditions (RH < 50) by Alexandris et al., (S9.11)
  }
  
  PET.Daily <- ET_truc.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Truc"
  PET_type <- "Reference Crop Evapotranspiration"
  
  if (solar == "data") {
    message1 <- "Solar radiation data has been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (humid == TRUE) {
    message4 <- "Adjustment for non-humid conditions has been applied to calculated Truc reference crop evapotranspiration"
  } else {
    message4 <- "No adjustment for non-humid conditions has been applied to calculated Truc reference crop evapotranspiration"
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message4)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message4=message4)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.HargreavesSamani <- function(data, constants, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Hargreaves-Samani
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  C_HS <- 0.00185 * (data$Tmax - data$Tmin)^2 - 0.0433 * (data$Tmax - data$Tmin) + 0.4023 # empirical coefficient by Hargreaves and Samani (1985) (S9.13)
  ET_HS.Daily <- 0.0135 * C_HS * R_a / constants$lambda * (data$Tmax - data$Tmin)^0.5 * (Ta + 17.8) # reference crop evapotranspiration by Hargreaves and Samani (1985) (S9.12)

  PET.Daily <- ET_HS.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Hargreaves-Samani"
  PET_type <- "Reference Crop Evapotranspiration"
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: reference crop")
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------
  
Evapotranspiration.ChapmanAustralian <- function(data, constants, Penpan, solar, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (Penpan == TRUE) { # Calculate Class-A pan evaporation using Penpan formula
    if (is.null(data$Tmax)|is.null(data$Tmin)) { 
      stop("fail to obtain data of 'Ta' (average daily temperature), 'vabar' (mean daily actual vapour pressure) and 'vas' (daily saturated vapour pressure)")
    }
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("fail to obtain data of 'vabar' (mean daily actual vapour pressure)")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("fail to obtain data of wind speed")
    }
    if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
      stop("fail to obtain data of 'Rs' (daily solar radiation)")
    } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
      stop("fail to obtain data of 'n' (daily sunshine hours)")
    } else if (solar == "cloud" & is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
      stop("fail to obtain data of 'cloud.cover.in.oktas' (daily sunshine hours)")
    } else if (solar == "monthly precipitation" & is.null(data$Cd)) { # for alternative calculation of cloudiness using monthly precipitation
      stop("fail to obtain data of 'monthly.precipitation.in.mm' (monthly precipitation)")
    } 
  }
  if (Penpan == FALSE & is.null(data$Epan)) { # for using Class-A pan evaporation data
    stop("fail to obtain 'Epan' (daily Class-A pan evaporation)")
  }
  # check user-input albedo
  if (Penpan == TRUE) {
    if (is.na(as.numeric(alpha))) {
      stop("Please use a numeric value for the alpha (albedo of evaporative surface)")
    }
    if (!is.na(as.numeric(alpha))) {
      if (as.numeric(alpha) < 0 | as.numeric(alpha) > 1) {
        stop("Please use a value between 0 and 1 for the alpha (albedo of evaporative surface)")
      }
    } 
  }

  # Calculations from data and constants for daily equivalent Penman-Monteith potential evaporation 
  A_p <- 0.17 + 0.011 * abs(constants$lat) # constant (S13.2)
  B_p <- 10 ^ (0.66 - 0.211 * abs(constants$lat)) # constants (S13.3)
  
  if (Penpan == TRUE) {
    # Calculating mean temperature 
    Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
    
    # Saturated vapour pressure
    vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
    vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
    vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
    
    # Vapour pressure
    vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
    
    # estimating class-A pan evaporation using Penpan model 
    P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
    delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
    gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
    d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
    delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
    w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
    N <- 24/pi * w_s # calculating daily values
    R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
    R_so <- (0.75 + (2*10^-5)*constants$Elev) * R_a # clear sky radiation (S3.4)
    
    if (solar == "data") {
      R_s <- data$Rs
    } else if (solar == "monthly precipitation") {
      # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
      R_s <- (0.85 - 0.047*data$Cd)*R_a 
    } else {
      # calculate R_s from sunshine hours - data or estimation using cloudness
      R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
    }
    
    # Wind speed
    if (is.null(data$u2)) {
      u2 <- data$uz * 4.87 / log(67.8*constants$z - 5.42) # Equation S5.20 for PET formulations other than Penman
    } else {
      u2 <- data$u2
    }
    
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For short grass
    P_rad <- 1.32 + 4 * 10^(-4) * abs(constants$lat) + 8 * 10^(-5) * (constants$lat)^2 # pan radiation factor (S6.6)
    f_dir <- -0.11 + 1.31 * R_s / R_a # fraction of R_S that os direct (S6.5)
    R_span <- (f_dir * P_rad + 1.42 * (1 - f_dir) + 0.42 * alpha) * R_s # total shortwave radiation received (S6.4)
    R_npan <-(1 - constants$alphaA) * R_span - R_nl # net radiation at the pan (S6.3)
    f_pan_u <-1.201 + 1.621 * u2 # (S6.2)
    
    Epan <- delta / (delta + constants$ap * gamma) * R_npan / constants$lambda + constants$ap * gamma / (delta + constants$ap * gamma) * f_pan_u * (vas - vabar) # Penpan estimation of Class-A pan evaporation (S6.1)
    ET_eqPM.Daily <- A_p * Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  } else if (Penpan == FALSE & is.null(data$Epan)) {
    stop("No data available for Class-A pan evaporation ")
  } else {
    ET_eqPM.Daily <- A_p * data$Epan + B_p # daily equivalent Penman-Monteith potential evaporation (mm.day^-1)
  }
 
  PET.Daily <- ET_eqPM.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (solar == "data") {
    message1 <- "Solar radiation data has been used for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
  }
  
  if (Penpan == TRUE) {
    message5 <- "Penpan formulation has been used to estimate Class-A pan evaporation for the calculation of equivalent Penman-Monteith potential evaporation"
  } else {
    message5 <- "Class-A pan evaporation has been used for the calculation of equivalent Penman-Monteith potential evaporation"
  }
  
  # Generate summary message for results
  PET_formulation <- "Chapman"
  PET_type <- "Equivalent Penmen-Monteith Reference Crop Evapotranspiration"
  if (Penpan == TRUE) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else {
    Surface <- paste("not specified, actual Class-A pan evaporation data is used")
  }
  
  message(PET_formulation, " ", PET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message5)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=message1, message5=message5)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

Evapotranspiration.JensenHaise <- function(data, constants, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Jensen-Haise
  # estimating evapotranspiration using Jensen-Haise
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  ET_JH.Daily <- R_a * Ta / (constants$lambda * 40) # Jensen-Haise daily evapotranspiration by Jensen and Haise  (1963) (mm.day^-1) (Oudin et al., 2005)
 
  PET.Daily <- ET_JH.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "Jensen-Haise"
  PET_type <- "Potential Evapotranspiration"
  
  message(PET_formulation, " ", PET_type)

  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type)
  class(results) <- funname
  
  return(results)
}

#-------------------------------------------------------------------------------------

Evapotranspiration.McGuinnessBordne <- function(data, constants, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("fail to obtain data of 'Ta' (average daily temperature)")
  }

  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for McGuinness-Bordne
  # estimating evapotranspiration using McGuinness-Bordne
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  ET_MB.Daily <- R_a * (Ta + 5) / (constants$lambda * 68) # McGuinness-Bordne daily evapotranspiration by McGuinness-Bordne  (1972) (mm.day^-1) (Oudin et al., 2005)
  
  PET.Daily <- ET_MB.Daily
  PET.Monthly <- aggregate(PET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  PET.Annual <- aggregate(PET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  PET.MonthlyAve <- PET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    PET.MonthlyAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    PET.AnnualAve[i] <- mean(PET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  PET_formulation <- "McGuinness-Bordne"
  PET_type <- "Potential Evapotranspiration"
  
  message(PET_formulation, " ", PET_type)
  
  results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type)
  class(results) <- funname
  
  return(results)
}
  #####################################################

  # Calculate radiation variables
Radiation <- function(data, constants, solar, Tdew) {
  class(data) <- funname

    # Check of specific data requirement
    if (is.null(data$Tmax)) {
      stop("fail to obtain data of 'Tmax' (daily maximum temperature)")
    }
    if (is.null(data$Tmin)) {
      stop("fail to obtain data of 'Tmin' (daily minimum temperature)")
    }
    if (Tdew == TRUE & is.null(data$Tdew)) {
      stop("fail to obtain data of 'Tdew' (daily dew point temperature")
    }
    if (Tdew == FALSE & (is.null(data$RHmax)|is.null(data$RHmin))) {
      stop("fail to obtain data of 'vabar' mean daily actual vapour pressure")
    }
    if (is.null(data$n)) {
      stop("fail to obtain data of 'n' (daily sunshine hours)")
    }
    if (solar == "data" | solar == "monthly precipitation") {
      stop("Only 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
    }  
  
    if (is.null(data$Precip)) {
      if ("PA" %in% names(constants) == FALSE) { 
        stop("fail to obtain data of 'PA' (annual average rainfall)")
      } # if annual average rainfall is not in data check the constants file
    } 
  
    # Convert daily data to required monthly data
    Tmax_Mo <- aggregate(data$Tmax, as.yearmon(data$Date.daily, "%m/%y"), FUN = max)
    Tmin_Mo <- aggregate(data$Tmin, as.yearmon(data$Date.daily, "%m/%y"), FUN = min)
    T_Mo <- (Tmax_Mo + Tmin_Mo) / 2
    
    # Dew point temperature 
    if (Tdew == TRUE) {
      Tdew_Mo <- aggregate(data$Tdew, as.yearmon(data$Date.daily, "%m/%y"), FUN = mean)
    } else {
      # Saturated vapour pressure
      vs_Tmax <- 0.6108 * exp(17.27 * data$Tmax / (data$Tmax + 237.3)) # Equation S2.5
      vs_Tmin <- 0.6108 * exp(17.27 * data$Tmin / (data$Tmin + 237.3)) # Equation S2.5
      vas <- (vs_Tmax + vs_Tmin)/2 # Equation S2.6
      
      # Vapour pressure
      vabar <- (vs_Tmin * data$RHmax/100 + vs_Tmax * data$RHmin/100)/2 # Equation S2.7
      
      vabar_Mo <- aggregate(vabar, as.yearmon(data$Date.daily, "%m/%y"), FUN = mean)
      Tdew_Mo <- (116.9 + 237.3*log(vabar_Mo)) / (16.78 - log(vabar_Mo))
    }
    
    # calculating ratio of daily sunshine hours to maximum daily sunshine hours, according to Morton's, averaged over all records
    
    # Slope of saturation vapour pressure curve [kPa/C]
    delta <- 4098 * (0.6108 * exp(17.27 * T_Mo / (T_Mo + 237.3)))/ (T_Mo + 237.3)^2 # Equation S2.4 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 13 in Allen et al, 1998. 
    
    # Mean daily maximum sunshine hours
    deltas <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar declination (rad) (S3.7)
    omegas <- acos(-tan(constants$lat_rad) * tan(deltas)) # sunset hour angle (rad) (S3.8)
    
    N <- 24/pi * omegas # calculating daily values
    
    S_daily <- data$n/N # daily ratios
    for (i in 1:length(S_daily)) {
      if (S_daily[i] > 1) {
        S_daily[i] <- 1
      }
    } # Criteria daily sunshine hours <= maximum daily sunshine hours
    
    S <- mean(S_daily) # according to Morton's, averaged over all records
    
    # Annual average rainfall
    if ("PA" %in% names(constants) == TRUE) {
      PA <- constants$PA
    } else {
      PA <- mean(aggregate(data$Precip, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum))
    }
    
                    
    # Update the following constant values for CRWE and CRLE  
    if (class(data) == "MortonCRLE" | class(data) =="MortonCRWE") {
      constants$epsilonMo <- 0.97 # (Morton, 1983)
      constants$fz <- 25.0 # Wm^-2.mbar^-1 for T >= 0 degree Celcius (Morton, 1983)
      constants$b0 <- 1.12 # (Morton, 1983)
      constants$b1 <- 13 # W.m^-2 (Morton, 1983)
      constants$b2 <- 1.12 # (Morton, 1983)
    }
    # calculate radiation
    ptops <- ((288 - 0.0065 * constants$Elev) / 288)^5.256 # ratio of atmospheric pressure to sea-level pressure (S21.3)
    alpha_zd <- 0.26 - 0.00012 * PA * sqrt(ptops) * (1 + abs(constants$lat/42) + (constants$lat/42)^2) # zenith value of the dry-season snow-free clear sky albedo 
    if (alpha_zd < 0.11) { 
      alpha_zd <- 0.11 
    } else {
      if (alpha_zd > 0.17) {
        alpha_zd <- 0.17
      } else {
        alpha_zd <- alpha_zd
      }
    } # constraint 0.11 <= alpha_zd <= 0.17 (S21.7)
    vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo / (Tdew_Mo + constants$betaMo)) # Saturation vapour pressure at dew point temperature (S21.8)
    v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo / (T_Mo + constants$betaMo)) # Saturation vapour pressure at air temperature (S21.10)
    deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo+constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.12)
    thetaMo <- (23.2 * sin((29.5 * data$i - 94) * pi/180)) * pi/180 # angle of extra-terrestrial global radiation (S21.14)
    Z_Mo <- acos(cos(constants$lat_rad - thetaMo)) # (S21.16)
    for (i in 1:length(Z_Mo)) {
      if (cos(Z_Mo[i]) < 0.001) {
        Z_Mo[i] <- acos(0.001)
      } 
    }# constraint cosZ_Mo >= 0.001 (S21.19)
    omegaMo <- acos(1 - cos(Z_Mo)/(cos(constants$lat_rad)*cos(thetaMo))) # (S21.20)
    cosz <- cos(Z_Mo) + (sin(omegaMo)/omegaMo - 1) * cos(constants$lat_rad) * cos(thetaMo) # (S21.23)
    etaMo <- 1 + 1/60*sin((29.5 * data$i - 106) * pi/180) # (S21.25)
    G_E <- 1354/(etaMo^2) * omegaMo/pi * cosz # Extra-terrestrial blobal radiation (S21.27)
    alpha_zz <- matrix(NA,length(v_Mo),1)
    alpha_zz[1:length(v_Mo)] <- alpha_zd # zenith value of snow-free clear sky albedo (S21.29)
    for (i in 1:length(v_Mo)) {
      if (alpha_zz[i] < 0.11) { 
        alpha_zz[i] <- 0.11 
      } else {
        if (alpha_zz[i] > 0.5 * (0.91 - vD_Mo/v_Mo)) {
          alpha_zz[i] <- 0.91 - vD_Mo/v_Mo
        } else {
        }
      } # constraint 0.11 <= alpha_zz <= 0.5 * (0.91 - vD_Mo/v_Mo) (S21.30)
    }
    
    # Update the alpha_zz value for CRWE and CRLE
    if (class(data) == "MortonCRLE" | class(data) =="MortonCRWE") {
      alpha_zz[1:length(v_Mo)] <- 0.05
    }
    
    c_0 <- as.vector(v_Mo - vD_Mo) # (S21.32)
    for (i in 1:length(c_0)) {
      if (c_0[i] < 0) { 
        c_0[i] <- 0
      } else {
        if (c_0[i] > 1) {
          c_0[i] <- 1
        } else {
          c_0[i] <- c_0[i]
        }
      } # constraint 0 <= c_0 <= 1 (S21.34) 
    }
    alpha_z <- alpha_zz + (1 - c_0^2) * (0.34 - alpha_zz) # zenith value of clear-sky albedo (S21.35)
    alpha_0 <- alpha_z * (exp(1.08) - ((2.16 * cos(Z_Mo))/pi + sin(Z_Mo)) * exp(0.012 * Z_Mo * 180/pi)) / (1.473 * (1 - sin(Z_Mo))) # clear-sky albedo (S21.37)
    W_Mo <- vD_Mo/(0.49 + T_Mo/129) # mm, precipirable water (S21.39)
    c_1 <- as.vector(21 - T_Mo) #(S21.41)
    for (i in 1:length(c_1)) {
      if (c_1[i] < 0) {
        c_1[i] <- 0
      } else {
        if (c_1[i] > 5) {
          c_1[i] <- 5
        } else {
          c_1[i] <- c_1[i]
        }
      } 
    }# constraint 0 <= c_1 <= 5 (S21.43)
    j_Mo <- (0.5 + 2.5 * (cosz)^2) * exp(c_1 * (ptops - 1)) # turbidity coefficient (S21.44)
    tauMo <- exp(-0.089 * (ptops * 1/cosz)^0.75 - 0.083 * (j_Mo/cosz)^0.9 - 0.029*(W_Mo/cosz)^0.6) # transmittancy of clear sky to direct bean radiation (S21.46)
    tauaMo <- as.vector(exp(-0.0415 * (j_Mo/cosz)^0.9 - (0.0029)^0.5 * (W_Mo/cosz)^0.3)) # proportion of tau_Mo that is the result of absorption
    for (i in 1:length(tauaMo)) {
      if (tauaMo[i] < exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 - 0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)) {
        tauaMo[i] <- exp(-0.0415 * (as.matrix(j_Mo/cosz)[i])^0.9 - 0.029 * (as.matrix(W_Mo/cosz)[i])^0.6)
      } else {
        tauaMo[i] <- tauaMo[i]
      }
    }# constraint tauaMo >= exp(-0.0415 * (j_Mo/cosz)^0.9 - 0.029 * (W_Mo/cosz)^0.6) (S21.51)
    G_0 <- G_E * tauMo * (1 + (1 - tauMo/tauaMo) * (1 + alpha_0 * tauMo)) # Wm^-2, clear-sky global radiation (S21.52)
    G_Mo <- S * G_0 + (0.08 + 0.30 * S) * (1 - S) * G_E # Wm^-2, incident global radiation (S21.54)
    alpha_Mo <- alpha_0 * (S + (1 - S) * (1 - Z_Mo/330 * 180/pi)) # average albedo (S21.56)
    c_2 <- as.vector(10 * (vD_Mo/v_Mo - S - 0.42)) # (S21.59)
    for (i in 1:length(c_2)) {
      if (c_2[i] < 0) {
        c_2[i] <- 0
      } else {
        if (c_2[i] > 1) {
          c_2[i] <- 1
        } else {
          c_2[i] <- c_2[i]
        }
      }
    }# constraint 0 <= c_2 <= 1 
    rouMo <- 0.18 * ((1 - c_2) * (1 - S)^2 + c_2 * (1 - S)^0.5) * 1/ptops # proportional increase in atmospheric radiation due to clouds (S21.60)
    B_Mo <- as.vector(constants$epsilonMo * constants$sigmaMo * (T_Mo +273)^4 * (1 - (0.71 + 0.007 * vD_Mo * ptops) * (1 + rouMo))) # Wm^-2, net longwave radiation loss for soil-plant surface at air temperature (S21.62)
    for (i in 1:length(B_Mo)) {
      if (B_Mo[i] < 0.05 * constants$epsilonMo * constants$sigmaMo * (T_Mo[i] + 274)^4) {
        B_Mo[i] <- 0.05 * constants$epsilonMo * constants$sigmaMo * (T_Mo[i] + 274)^4
      } else {
        B_Mo[i] <- B_Mo[i]
      } 
    }# constraint B_Mo < 0.05 * constants$epsilonMo * sigmaMo * (T_Mo + 274)^4 (S21.64)
  
    # Generate summary message for results
    if (solar == "data") {
      message1 <- "Sunshine hour data has been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
      message1 <- "Cloudiness data has been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
      message1 <- "Monthly precipitation data has been used for calculating incoming solar radiation"
    }
  
    if (Tdew == TRUE) {
      message6 <- "Data of dew point temperature has been used"
    } else {
      message6 <- "Data of average vapour pressure has been used to estimate dew point pressure"
    }
  
    variables <- list(Tmax_Mo=Tmax_Mo, Tmin_Mo=Tmin_Mo, T_Mo=T_Mo, Tdew_Mo=Tdew_Mo, S=S, ptops=ptops, vD_Mo=vD_Mo, v_Mo=v_Mo, deltaMo=deltaMo, G_E=G_E, G_Mo=G_Mo, alpha_Mo=alpha_Mo, B_Mo=B_Mo, message1=message1, message6=message6)
    return(variables)
  }
  
  #-------------------------------------------------------------------------------------
  Evapotranspiration.MortonCRAE <- function(data, constants, est, solar, Tdew, ...)  {
    
    variables <- Radiation(data, constants, solar, Tdew)
    
    # Morton's CRAE procedure
    
    R_T <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
    R_TC <- as.vector(R_T) # Wm^-2, (S21.68)
    for (i in 1:length(R_TC)) {
      if (R_TC[i] < 0) {
        R_TC[i] <- 0
      } else {
        R_TC[i] <- R_TC[i]
      } # constraint R_TC >= 0 
    }
    xiMo <- 1/(0.28 * (1 + variables$vD_Mo/variables$v_Mo) + R_TC * variables$deltaMo / (variables$ptops * constants$gammaps * (1/variables$ptops)^0.5 * constants$b0 * constants$fz * (variables$v_Mo - variables$vD_Mo))) # a dimensionless stability factor (S21.69)
    for (i in 1:length(xiMo)) {
      if (xiMo[i] < 1) {
        xiMo[i] <- 1
      } else {
        xiMo[i] <- xiMo[i]
      } # constraint xiMo >= 1 
    }
    f_T <- (1/variables$ptops)^0.5 * constants$fz / xiMo # vapour transfer coefficient (S21.71)
    lambdaMo1 <- constants$gammaps * variables$ptops + 4 * constants$epsilonMo * constants$sigmaMo * (variables$T_Mo + 274)^3 / f_T # heat transfer coefficient (S21.73)
    # Iteration for equilibrium temperature T_p
    T_p <- variables$T_Mo
    for (i in 1:99999) {
      v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
      delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
      delta_T_p <- (R_T/f_T + variables$vD_Mo - v_p + lambdaMo1 * (variables$T_Mo - T_p)) / (delta_p + lambdaMo1) # change in T_p (S21.75)
      T_p <- T_p + delta_T_p # T_p for next iteration (S21.76)
      if (abs(max(na.omit(delta_T_p))) < 0.01) break
    }
    v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
    delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
    
    # Apply Morton Potential Point Evaporation
    E_TP.temp <- R_T - lambdaMo1 * f_T * (T_p - variables$T_Mo) # Wm^-2, potential evapotranspiration (S21.79)
    # Apply Morton Potential Wet Evaporation
    R_TP <- E_TP.temp + variables$ptops * constants$gammaps * f_T * (T_p - variables$T_Mo) # Wm^-2, net radiation at the soil-plant surface for equilibrium temperature (S21.81)
    E_TW.temp <- constants$b1 + constants$b2 * R_TP / (1 + variables$ptops * constants$gammaps / delta_p) # Wm^-2, wet-environment areal evapotranspiration (S21.84)
    # Apply Morton Potential Areal Evaporation
    E_T_Mo.temp <- 2 * E_TW.temp - E_TP.temp # Wm^-2, actual areal evapotranspiration (S21.86)
    
    # Convert evaporation in power unit of W.m^-2 to evaporation units of mm.day^-1
    E_TP.temp <- 1/(constants$lambdaMo) * E_TP.temp # mm.day^-1 (S21.88)
    E_TW.temp <- 1/(constants$lambdaMo) * E_TW.temp # mm.day^-1 (S21.89)
    E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp # mm.day^-1 (S21.90)
    
    # Calculate monthly evaporation in mm.month^-1
    E_TP <- E_TP.temp * data$ndays
    E_TW <- E_TW.temp * data$ndays
    E_T_Mo <- E_T_Mo.temp * data$ndays
    
    if (est == "potential ET") {
      ET_Mo.Monthly <- E_TP
      ET_Mo.Average <- E_TP.temp
      PET_type <- "Potential Evapotranspiration"
    } else if (est == "wet areal ET") {
      ET_Mo.Monthly <- E_TW
      ET_Mo.Average <- E_TW.temp
      PET_type <- "Wet-environment Areal Evapotranspiration"
    } else if (est == "actual areal ET") {
      ET_Mo.Monthly <- E_T_Mo
      ET_Mo.Average <- E_T_Mo.temp
      PET_type <- "Actual Areal Evapotranspiration"
    }
    
    PET.Daily <- NULL
    PET.Monthly <- ET_Mo.Monthly
    PET.Annual <- aggregate(PET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
    
    PET.MonthlyAve <- PET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)){
      i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
      PET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon== mon])
    }
    for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)){
      i = year - min(as.POSIXlt(data$Date.monthly)$year) + 1
      PET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year== year])
    }
    
    # Generate summary message for results
    PET_formulation <- "Morton CRAE"
    
    message(PET_formulation, " ", PET_type)
    message(variables$message1)
    message(variables$message6)

    results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=variables$message1, message6=variables$message6)
    class(results) <- funname
    
    return(results)
    # End of Morton's CRAE procedure
  }
  
  
  #-----------------------------------------------------------------------------------
  Evapotranspiration.MortonCRWE <- function(data, constants, est, solar, Tdew, ...) {

    constants$epsilonMo <- 0.97 # (Morton, 1983)
    constants$fz <- 25.0 # Wm^-2.mbar^-1 for T >= 0 degree Celcius (Morton, 1983)
    constants$b0 <- 1.12 # (Morton, 1983)
    constants$b1 <- 13 # W.m^-2 (Morton, 1983)
    constants$b2 <- 1.12 # (Morton, 1983)
    variables <- Radiation(data, constants, solar, Tdew)
    
    # Morton's CRWE procedure
    
    alpha_zz <- 0.05
    
    R_W <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
    R_TC <- as.vector(R_W) # Wm^-2, (S21.68)
    for (i in 1:length(R_TC)) {
      if (R_TC[i] < 0) {
        R_TC[i] <- 0
      } else {
        R_TC[i] <- R_TC[i]
      } # constraint R_TC >= 0 
    }
    xiMo <- 1/(0.28 * (1 + variables$vD_Mo/variables$v_Mo) + R_TC * variables$deltaMo / (variables$ptops * constants$gammaps * (1/variables$ptops)^0.5 * constants$b0 * constants$fz * (variables$v_Mo - variables$vD_Mo))) # a dimensionless stability factor (S21.69)
    for (i in 1:length(xiMo)) {
      if (xiMo[i] < 1) {
        xiMo[i] <- 1
      } else {
        xiMo[i] <- xiMo[i]
      } # constraint xiMo >= 1 
    }
    f_T <- (1/variables$ptops)^0.5 * constants$fz / xiMo # vapour transfer coefficient (S21.71)
    lambdaMo1 <- constants$gammaps * variables$ptops + 4 * constants$epsilonMo * constants$sigmaMo * (variables$T_Mo + 274)^3 / f_T # heat transfer coefficient (S21.73)
    # Iteration for equilibrium temperature T_p
    T_p <- variables$T_Mo
    for (i in 1:99999) {
      v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
      delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
      delta_T_p <- (R_W/f_T + variables$vD_Mo - v_p + lambdaMo1 * (variables$T_Mo - T_p)) / (delta_p + lambdaMo1) # change in T_p (S21.75)
      T_p <- T_p + delta_T_p # T_p for next iteration (S21.76)
      if (abs(max(na.omit(delta_T_p))) < 0.01) break
    }
    v_p <- 6.11 * exp((constants$alphaMo * T_p)/(T_p + constants$betaMo)) # mbar, saturation vapour pressure at equilibrium temperature (S21.77)
    delta_p <- constants$alphaMo * constants$betaMo * v_p/((T_p + constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.78)
    
    # Apply Morton Potential Point Evaporation
    E_P.temp <- R_W - lambdaMo1 * f_T * (T_p - variables$T_Mo) # Wm^-2, potential evapotranspiration (S21.79)
    # Apply Morton Potential Wet Evaporation
    R_P <- E_P.temp + variables$ptops * constants$gammaps * f_T * (T_p - variables$T_Mo) # Wm^-2, net radiation at the water surface for equilibrium temperature (S21.81)
    E_W.temp <- constants$b1 + constants$b2 * R_P / (1 + variables$ptops * constants$gammaps / delta_p) # Wm^-2, wet-environment areal evapotranspiration (S21.84)
    # Apply Morton Potential Areal Evaporation
    E_T_Mo.temp <- 2 * E_W.temp - E_P.temp # Wm^-2, actual areal evapotranspiration (S21.86)
    
    # Convert evaporation in power unit of W.m^-2 to evaporation units of mm.day^-1
    E_P.temp <- 1/(constants$lambdaMo) * E_P.temp # mm.day^-1 (S21.88)
    E_W.temp <- 1/(constants$lambdaMo) * E_W.temp # mm.day^-1 (S21.89)
    E_T_Mo.temp <- 1/(constants$lambdaMo) * E_T_Mo.temp # mm.day^-1 (S21.90)
    
    # Calculate monthly evaporation in mm.month^-1
    E_P <- E_P.temp * data$ndays
    E_W <- E_W.temp * data$ndays
    E_T_Mo <- E_T_Mo.temp * data$ndays
    
    if (est == "potential") {
      ET_Mo.Monthly <- E_P
      ET_Mo.Average <- E_P.temp
      PET_type <- "Potential Evaporation"
    } else if (est == "shallow lake") {
      ET_Mo.Monthly <- E_W
      ET_Mo.Average <- E_W.temp
      PET_type <- "Shallow Lake Evaporation"
    }

    PET.Daily <- NULL
    PET.Monthly <- ET_Mo.Monthly
    PET.Annual <- aggregate(PET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
    
    PET.MonthlyAve <- PET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)){
      i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
      PET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon== mon])
    }
    for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)){
      i = year - min(as.POSIXlt(data$Date.monthly)$year) + 1
      PET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year== year])
    }
    
    # Generate summary message for results
    PET_formulation <- "Morton CRWE"
    
    message(PET_formulation, " ", PET_type)
    message(variables$message1)
    message(variables$message6)

    results <- list(PET.Daily=PET.Daily, PET.Monthly=PET.Monthly, PET.Annual=PET.Annual, PET.MonthlyAve=PET.MonthlyAve, PET.AnnualAve=PET.AnnualAve, PET_formulation=PET_formulation, PET_type=PET_type, message1=variables$message1, message6=variables$message6)
    class(results) <- funname
    
    return(results)
  }
  # End of Morton's CRWE procedure
  
  
  #-------------------------------------------------------------------------------------

