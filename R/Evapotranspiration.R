
ET <- function(data, ...) UseMethod("ET")

  #-------------------------------------------------------------------------------------
  
ET.Penman <- function(data, constants, solar, wind, windfunction_ver, alpha = 0.08, z0 = 0.001, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
    }
    if (is.null(data$uz) & is.null(data$u2)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }

  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
   
   ET.Daily <- Epenman.Daily
   ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
   ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
   ET.MonthlyAve <- ET.AnnualAve <- NULL
   for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
     i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
     ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
   }
   for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
     i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
     ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
   }

   # Generate summary message for results  
   ET_formulation <- "Penman" 
   if (wind == "no") {
     ET_type <- "Open-water Evaporation"
     Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
   } else {
     if (alpha != 0.08) {
       ET_type <- "Potential ET"
       Surface <- paste("user-defined, albedo =", alpha, "; roughness height =", z0, "m")
     } else if (alpha == 0.08) {
       ET_type <- "Open-water Evaporation"
       Surface <- paste("water, albedo =", alpha, "; roughness height =", z0, "m")
     }
   }
   
   if (solar == "data") {
     message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"    
   } else if (solar == "sunshine hours") {
     message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
   } else if (solar == "cloud") {
     message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
   } else {
     message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
   }
   
   if (wind == "yes") {
     if (windfunction_ver == "1948") {
       message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
     } else if (windfunction_ver == "1956") {
       message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
     } 
   } else {
     message2 <- "Alternative calculation for Penman evaporation without wind data have been performed"
   }
  
   message(ET_formulation, " ", ET_type)
   message("Evaporative surface: ", Surface)
   message(message1)
   message(message2)
  
   results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
   class(results) <- funname
   return(results)
}

  #-------------------------------------------------------------------------------------

ET.PenmanMonteith <- function(data, constants, solar, wind, crop, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
  stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
    ET_RC.Daily <- matrix(NA,length(data$date.Daily),1)
    ET_RC.Monthly <- 0.038 * R_s.Monthly * sqrt(Ta.Monthly + 9.5) - 2.4 * (R_s.Monthly/R_a.Monthly)^2 + 0.075 * (Ta.Monthly + 20) * (1 - RHmean.Monthly/100) # Reference crop evapotranspiration without wind data by Valiantzas (2006) (S5.21)
  }
  
  ET.Daily <- ET_RC.Daily
  if (is.na(mean(ET_RC.Daily))) {
    ET_RC.Daily <- data$Tmax
    for (cont in 1:length(data$i)) {
      ET_RC.Daily[(((as.numeric(as.yearmon(time(ET_RC.Daily))))-floor(as.numeric(as.yearmon(time(ET_RC.Daily)))))*12+1)==data$i[cont]] <- ET_RC.Monthly[cont]
    }
    ET.Daily <- ET_RC.Daily
    ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
    ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
  } else {
    ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
    ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  }
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (wind == "no") {
    ET_formulation <- "Penman-Monteith (without wind data)"
    ET_type <- "Reference Crop ET"
    Surface <- paste("short grass, albedo =", alpha, "; roughness height =", z0, "m")
  } else {
    if (crop == "short") {
      ET_formulation <- "Penman-Monteith FAO56"
      ET_type <- "Reference Crop ET"
      Surface <- paste("FAO-56 hypothetical short grass, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    } else {
      ET_formulation <- "Penman-Monteith ASCE-EWRI Standardised"
      ET_type <- "Reference Crop ET"
      Surface <- paste("ASCE-EWRI hypothetical tall grass, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, " m; roughness height =", z0, "m")
    }
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    message2 <- "Wind data have been used for calculating the reference crop evapotranspiration"
  } else {
    message2 <- "Alternative calculation for reference crop evapotranspiration without wind data have been performed"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  class(results) <- funname
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.MattShuttleworth <- function(data, constants, solar, alpha, r_s, CH, ...) { 
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
  r_clim <- 86400 * constants$Roua * constants$Ca * (vas - vabar) / (delta * R_ng) # clinmatological resistance (s*m^-1) (S5.34)
  r_clim[r_clim == 0] <- 0.1 # correction for r_clim = 0
  u2[u2 == 0] <- 0.1 # correction for u2 = 0
  VPD50toVPD2 <- (302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) + 1/r_clim * ((302 * (delta + gamma) + 70 * gamma * u2) / (208 * (delta + gamma) + 70 * gamma * u2) * (208 / u2) - (302 / u2)) # ratio of vapour pressure deficits at 50m to vapour pressure deficits at 2m heights (S5.35)
  r_c50 <- 1 / ((0.41)^2) * log((50 - 0.67 * CH) / (0.123 * CH)) * log((50 - 0.67 * CH) / (0.0123 * CH)) * log((2 - 0.08) / 0.0148) / log((50 - 0.08) / 0.0148) # aerodynamic coefficient (s*m^-1) (S5.36)
  
  E_Tc.Daily <- 1/constants$lambda * (delta * R_ng + (constants$Roua * constants$Ca * u2 * (vas - vabar)) / r_c50 * VPD50toVPD2) / (delta + gamma * (1 + r_s * u2 / r_c50)) # well-watered crop evapotranspiration in a semi-arid and windy location (S5.37)
  
  ET.Daily <- E_Tc.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Matt-Shuttleworth"
  ET_type <- "Reference Crop ET"
  Surface <- paste("user-defined, albedo =", alpha, "; surface resisitance =", r_s, "sm^-1; crop height =", CH, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.PriestleyTaylor <- function(data, constants, solar, alpha, ...) {  
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
  ET.Daily <- E_PT.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Priestley-Taylor"
  ET_type <- "Potential ET"
  if (alpha != 0.08) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else if (alpha == 0.08) {
    Surface <- paste("water, albedo =", alpha)
  }
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.Penpan <- function(data, constants, solar, alpha, overest, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
    
  ET.Daily <- Epenpan.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Penpan"
  ET_type <- "Class-A Pan Evaporation"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.BrutsaertStrickler <- function(data, constants, solar, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
  ET.Daily <- ET_BS_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Brutsaert-Strickler"
  ET_type <- "Actual Areal ET"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.GrangerGray <- function(data, constants, solar, windfunction_ver, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (is.null(data$u2) & is.null(data$uz)) {
    stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
  ET.Daily <- ET_GG_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Granger-Gray"
  ET_type <- "Actual Areal ET"
  Surface <- paste("user-defined, albedo =", alpha)
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (windfunction_ver == "1948") {
    message2 <- "Wind data have been used for the calculation of the drying power of air, using Penman 1948 wind function."
  } else if (windfunction_ver == "1956") {
    message2 <- "Wind data have been used for the calculation of the drying power of air, using Penman 1956 wind function."
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)

  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.SzilagyiJozsa <- function(data, constants, solar, wind, windfunction_ver, alpha, z0, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (is.null(data$RHmax)|is.null(data$RHmin)) {
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  }
  if (wind == "yes") { # wind data is required
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
  }
  
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
    R_nl <- constants$sigma * (0.34 - 0.14 * sqrt(vabar)) * ((data$Tmax+273.2)^4 + (data$Tmin+273.2)^4)/2  * (1.35 * R_s / R_so - 0.35) # estimated net outgoing longwave radiation (S3.3)
    # For vegetated surface
    R_nsg <- (1 - alpha) * R_s # net incoming shortwave radiation (S3.2)
    R_ng = R_nsg - R_nl # net radiation (S3.1)
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
  
  ET.Daily <- E_SJ_Act.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Szilagyi-Jozsa"
  ET_type <- "Actual ET"
  Surface <- paste("user-defined, albedo =", alpha, "; roughness height", z0, "m")
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (wind == "yes") {
    if (windfunction_ver == "1948") {
      message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1948 wind function has been used."
    } else if (windfunction_ver == "1956") {
      message2 <- "Wind data have been used for calculating the Penman evaporation. Penman 1956 wind function has been used."
    } 
  } else {
    message2 <- "Alternative calculation for Penman evaporation without wind data have been performed"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message2)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message2=message2)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.Makkink <- function(data, constants, solar, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
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
  
  ET.Daily <- Emakkink.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Makkink"
  ET_type <- "Potential ET"

  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: open-water")
  message(message1)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.BlaneyCriddle <- function(data, constants, solar, height, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "sunshine hours" & is.null(data$n)) { # sunshine hour data is required
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud") {
    if (is.null(data$n)) { # for alternative calculation of sunshine hours using cloud cover
      stop("Required data missing for 'Cd.daily'")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
    if (is.null(data$RHmin)) { 
      stop("Required data missing for 'RHmin.daily'")
    } 
  } 
  if (solar == "data" | solar == "monthly precipitation") {
    stop("Only 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
  }
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  # Calculations from data and constants for Blaney and Criddle
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
  # chech if data from first/last years is incomplete, and adjust N.annual values for incomplete years
  if(data$J[1]!=1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4)==FALSE) { # first year a normal year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:365) - 1.39))))
  }
  if(data$J[1]!=1 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[1]/4)==TRUE) { # first year a leap year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[1]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:366) - 1.39))))
  }
  if (data$J[length(data$J)]!=365 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4)==FALSE) { # last year a normal year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/365 * c(1:365) - 1.39))))
  }
  if (data$J[length(data$J)]!=366 & is.integer(floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]/4)==TRUE) { # first year a leap year
    N.annual[floor(as.numeric(as.yearmon(data$Date.daily)))==floor(as.numeric(as.yearmon(data$Date.daily)))[length(data$J)]] <- 
      sum(24/pi * acos(-tan(constants$lat_rad) * tan(0.409 * sin(2*pi/366 * c(1:366) - 1.39))))
  }
  p_y <- 100 * data$n/N.annual # percentage of actual daytime hours for the day comparing to the annual sum of maximum sunshine hours

  
  ET_BC.Daily <- (0.0043 * data$RHmin - data$n/N - 1.41) + bvar * p_y * (0.46 * Ta +8.13) # Blaney-Criddle Reference Crop evapotranspiration (mm.day^-1) (S9.7)
  
  ET.Daily <- ET_BC.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  if (height == TRUE) {
    ET_BC.Daily = ET_BC.Daily * (1 + 0.1 * constants$Elev/1000) # with adjustment for site elevation by Allen and Pruitt (1986) (S9.9) 
  }
  
  # Generate summary message for results
  ET_formulation <- "Blaney-Criddle"
  ET_type <- "Reference Crop ET"

  if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (height == TRUE) {
    message3 <- "Height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
  } else {
    message3 <- "No height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message3)
    
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message3=message3)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.Turc <- function(data, constants, solar, humid, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  }

  if (humid == TRUE & (is.null(data$RHmax)|is.null(data$RHmin))) { # for adjustment for non-humid conditions
    stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
  } 
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  # Calculations from data and constants for Turc
  
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
  
  ET_Turc.Daily <- 0.013 * (23.88 * R_s + 50) * Ta / (Ta + 15) # reference crop evapotranspiration by Turc (1961) (S9.10)
  
  if (humid == TRUE) {
    # mean relative humidity
    RHmean <- (data$RHmax + data$RHmin) / 2 
    
    ET_Turc.Daily[RHmean < 50] <- 0.013 * (23.88 * R_s + 50) * Ta[RHmean < 50] / (Ta[RHmean < 50] + 15) * (1 + (50 - RHmean[RHmean < 50]) / 70) # Turc reference crop evapotranspiration adjusted for non-humid conditions (RH < 50) by Alexandris et al., (S9.11)
  }
  
  ET.Daily <- ET_Turc.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Turc"
  ET_type <- "Reference Crop ET"
  
  if (solar == "data") {
    message1 <- "Solar radiation data have been used directly for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (humid == TRUE) {
    message4 <- "Adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration"
  } else {
    message4 <- "No adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration"
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  message(message1)
  message(message4)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message4=message4)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.HargreavesSamani <- function(data, constants, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
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

  ET.Daily <- ET_HS.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Hargreaves-Samani"
  ET_type <- "Reference Crop ET"
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: reference crop")
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------
  
ET.ChapmanAustralian <- function(data, constants, Penpan, solar, alpha, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (Penpan == TRUE) { # Calculate Class-A pan evaporation using Penpan formula
    if (is.null(data$Tmax)|is.null(data$Tmin)) { 
      stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
    }
    if (is.null(data$RHmax)|is.null(data$RHmin)) {
      stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
    }
    if (is.null(data$u2) & is.null(data$uz)) {
      stop("Required data missing for 'uz.subdaily' or 'u2.subdaily'")
    }
    if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
      stop("Required data missing for 'Rs.daily'")
    } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
      stop("Required data missing for 'n.daily'")
    } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
      stop("Required data missing for 'Cd.daily'")
    } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
      stop("Required data missing for 'Precip.daily'")
    } 
  }
  if (Penpan == FALSE & is.null(data$Epan)) { # for using Class-A pan evaporation data
    stop("Required data missing for 'Epan.daily'")
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
 
  ET.Daily <- ET_eqPM.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  if (solar == "data") {
    message1 <- "Solar radiation data have been used for calculating evapotranspiration"
  } else if (solar == "sunshine hours") {
    message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
  } else if (solar == "cloud") {
    message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
  } else {
    message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
  }
  
  if (Penpan == TRUE) {
    message5 <- "Penpan formulation has been used to estimate Class-A pan evaporation for the calculation of equivalent Penman-Monteith potential evaporation"
  } else {
    message5 <- "Class-A pan evaporation has been used for the calculation of equivalent Penman-Monteith potential evaporation"
  }
  
  # Generate summary message for results
  ET_formulation <- "Chapman"
  ET_type <- "Equivalent Penmen-Monteith Reference Crop ET"
  if (Penpan == TRUE) {
    Surface <- paste("user-defined, albedo =", alpha)
  } else {
    Surface <- paste("not specified, actual Class-A pan evaporation data is used")
  }
  
  message(ET_formulation, " ", ET_type)
  message("Evaporative surface: ", Surface)
  message(message1)
  message(message5)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=message1, message5=message5)
  class(results) <- funname
  
  return(results)
}

  #-------------------------------------------------------------------------------------

ET.JensenHaise <- function(data, constants,solar,...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  if (solar == "data" & is.null(data$Rs)) { # solar radiation data is required
    stop("Required data missing for 'Rs.daily'")
  } else if (solar == "sunshine hours" & is.null(data$n)) { # for alternative calculation of solar radiation with sunshine hour
    stop("Required data missing for 'n.daily'")
  } else if (solar == "cloud" & is.null(data$Cd)) { # for alternative calculation of sunshine hours using cloud cover
    stop("Required data missing for 'Cd.daily'")
  } else if (solar == "monthly precipitation" & is.null(data$Precip)) { # for alternative calculation of cloudiness using monthly precipitation
    stop("Required data missing for 'Precip.daily'")
  } 
  
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 
  
  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)
  
  if (solar == "data") {
    R_s <- data$Rs
  } else if (solar == "monthly precipitation") {
    # calculate R_s from cloudness estimated from monthly precipitation (#S3.14)
    R_s <- (0.85 - 0.047*data$Cd)*R_a 
  } else {
    # calculate R_s from sunshine hours - data or estimation using cloudness
    R_s <- (constants$as + constants$bs * (data$n/N))*R_a # estimated incoming solar radiation (S3.9)
  }
  
  # estimating evapotranspiration using Jensen-Haise

  ET_JH.Daily <- (0.014 * (1.8*Ta+32) - 0.37) * R_s / constants$lambda # Jensen-Haise daily evapotranspiration by Jensen and Haise  (1963), *0.0394 for MJ.m^-2.day^-1 to in.day^-1, *33.8 for degree C to degree F 
 
  ET.Daily <- ET_JH.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "Jensen-Haise"
  ET_type <- "Potential ET"
  
  message(ET_formulation, " ", ET_type)

  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  class(results) <- funname
  
  return(results)
}

#-------------------------------------------------------------------------------------

ET.McGuinnessBordne <- function(data, constants, ...) {
  class(data) <- funname
  
  # Check of specific data requirement
  if (is.null(data$Tmax)|is.null(data$Tmin)) { 
    stop("Required data missing for 'Tmax.daily' and 'Tmin.daily', or 'Temp.subdaily'")
  }
  
  # Calculating mean temperature 
  Ta <- (data$Tmax + data$Tmin) / 2   # Equation S2.1 in Tom McMahon's HESS 2013 paper, which in turn was based on Equation 9 in Allen et al, 1998. 

  P <- 101.3 * ((293 - 0.0065 * constants$Elev) / 293)^5.26 # atmospheric pressure (S2.10)
  delta <- 4098 * (0.6108 * exp((17.27 * Ta)/(Ta+237.3))) / ((Ta + 237.3)^2) # slope of vapour pressure curve (S2.4)
  gamma <- 0.00163 * P / constants$lambda # psychrometric constant (S2.9)
  d_r2 <- 1 + 0.033*cos(2*pi/365 * data$J) # dr is the inverse relative distance Earth-Sun (S3.6)
  delta2 <- 0.409 * sin(2*pi/365 * data$J - 1.39) # solar dedication (S3.7)
  w_s <- acos(-tan(constants$lat_rad) * tan(delta2))  # sunset hour angle (S3.8)
  N <- 24/pi * w_s # calculating daily values
  R_a <- (1440/pi) * d_r2 * constants$Gsc * (w_s * sin(constants$lat_rad) * sin(delta2) + cos(constants$lat_rad) * cos(delta2) * sin(w_s)) # extraterristrial radiation (S3.5)

  # estimating evapotranspiration using McGuinness-Bordne

  ET_MB.Daily <- R_a * (Ta + 5)/ (constants$lambda*68)  # McGuinness-Bordne daily evapotranspiration by McGuinness-Bordne  (1972) (mm.day^-1) (Oudin et al., 2005
  
  ET.Daily <- ET_MB.Daily
  ET.Monthly <- aggregate(ET.Daily, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
  ET.Annual <- aggregate(ET.Daily, floor(as.numeric(as.yearmon(data$Date.daily, "%m/%y"))), FUN = sum)
  
  ET.MonthlyAve <- ET.AnnualAve <- NULL
  for (mon in min(as.POSIXlt(data$Date.daily)$mon):max(as.POSIXlt(data$Date.daily)$mon)){
    i = mon - min(as.POSIXlt(data$Date.daily)$mon) + 1
    ET.MonthlyAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$mon== mon])
  }
  for (year in min(as.POSIXlt(data$Date.daily)$year):max(as.POSIXlt(data$Date.daily)$year)){
    i = year - min(as.POSIXlt(data$Date.daily)$year) + 1
    ET.AnnualAve[i] <- mean(ET.Daily[as.POSIXlt(data$Date.daily)$year== year])
  }
  
  # Generate summary message for results
  ET_formulation <- "McGuinness-Bordne"
  ET_type <- "Potential ET"
  
  message(ET_formulation, " ", ET_type)
  
  results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type)
  class(results) <- funname
  
  return(results)
}
  #####################################################

  # Calculate radiation variables
Radiation <- function(data, constants, solar, Tdew, alpha=NULL) {
  class(data) <- funname

    # Check of specific data requirement
    if (is.null(data$Tmax)) {
      stop("Required data missing for 'Tmax.daily' or 'Temp.subdaily'")
    }
    if (is.null(data$Tmin)) {
      stop("Required data missing for 'Tmin.daily' or 'Temp.subdaily'")
    }
    if (Tdew == TRUE & is.null(data$Tdew)) {
      stop("Required data missing for 'Tdew.subdaily'")
    }
    if (Tdew == FALSE & (is.null(data$RHmax)|is.null(data$RHmin))) {
      stop("Required data missing for 'RHmax.daily' and 'RHmin.daily', or 'RH.subdaily'")
    }
    if (is.null(data$n)) {
      stop("Required data missing for 'n.daily'")
    }
    if (solar == "monthly precipitation") {
      stop("Only 'data', 'sunshine hours' and 'cloud' are accepted because estimations of sunshine hours is required")
    }  
  
    if (is.null(data$Precip)) {
      if ("PA" %in% names(constants) == FALSE) { 
        stop("Required data missing for 'Precip.daily' or required constant missing for 'PA'")
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
    
  if (solar == "sunshine hours") {
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
  } else if (solar == "data") {
    alpha_Mo = alpha
    vD_Mo <- 6.11 * exp(constants$alphaMo * Tdew_Mo / (Tdew_Mo + constants$betaMo)) # Saturation vapour pressure at dew point temperature (S21.8)
    v_Mo <- 6.11 * exp(constants$alphaMo * T_Mo / (T_Mo + constants$betaMo)) # Saturation vapour pressure at air temperature (S21.10)
    ptops <- ((288 - 0.0065 * constants$Elev) / 288)^5.256 # ratio of atmospheric pressure to sea-level pressure (S21.3)
    deltaMo <- constants$alphaMo * constants$betaMo * v_Mo/((T_Mo+constants$betaMo)^2) # mbar slope of vapour pressure curve (S21.12)
    G_E = NULL
    G_Mo = NULL
    S = NULL
    B_Mo = NULL
  }
    
    # Generate summary message for results
    if (solar == "sunshine hours") {
      message1 <- "Sunshine hour data have been used for calculating incoming solar radiation"
    } else if (solar == "cloud") {
      message1 <- "Cloudiness data have been used for calculating sunshine hour and thus incoming solar radiation"
    } else {
      message1 <- "Monthly precipitation data have been used for calculating incoming solar radiation"
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
  ET.MortonCRAE <- function(data, constants, est, solar, Tdew, alpha=NULL, ...)  {
    
    variables <- Radiation(data, constants, solar, Tdew)
    
    # Morton's CRAE procedure
    if (solar == "sunshine hours") {
      R_T <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
    } else if (solar == "data") {
      variables$alpha_Mo <- alpha
      Rs_Mo <- aggregate(data$Rs, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
      R_T <- (1 - variables$alpha_Mo) * Rs_Mo
    }

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
      ET_type <- "Potential ET"
    } else if (est == "wet areal ET") {
      ET_Mo.Monthly <- E_TW
      ET_Mo.Average <- E_TW.temp
      ET_type <- "Wet-environment Areal ET"
    } else if (est == "actual areal ET") {
      ET_Mo.Monthly <- E_T_Mo
      ET_Mo.Average <- E_T_Mo.temp
      ET_type <- "Actual Areal ET"
    }
    
    ET.Daily <- NULL
    ET.Monthly <- ET_Mo.Monthly
    ET.Annual <- aggregate(ET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
    
    ET.MonthlyAve <- ET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)){
      i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
      ET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon== mon])
    }
    for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)){
      i = year - min(as.POSIXlt(data$Date.monthly)$year) + 1
      ET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year== year])
    }
    
    # Generate summary message for results
    ET_formulation <- "Morton CRAE"
    
    message(ET_formulation, " ", ET_type)
    message(variables$message1)
    message(variables$message6)

    results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=variables$message1, message6=variables$message6)
    class(results) <- funname
    
    return(results)
    # End of Morton's CRAE procedure
  }
  
  
  #-----------------------------------------------------------------------------------
  ET.MortonCRWE <- function(data, constants, est, solar, Tdew, alpha=NULL, ...) {

    constants$epsilonMo <- 0.97 # (Morton, 1983)
    constants$fz <- 25.0 # Wm^-2.mbar^-1 for T >= 0 degree Celcius (Morton, 1983)
    constants$b0 <- 1.12 # (Morton, 1983)
    constants$b1 <- 13 # W.m^-2 (Morton, 1983)
    constants$b2 <- 1.12 # (Morton, 1983)
    variables <- Radiation(data, constants, solar, Tdew)
    
    # Morton's CRWE procedure
    
    alpha_zz <- 0.05
    
    # Morton's CRAE procedure
    if (solar == "sunshine hours") {
      R_W <- (1 - variables$alpha_Mo) * variables$G_Mo - variables$B_Mo # Wm^-2, net radiation at soil-plant surface at air temperature (S21.66)
    } else if (solar == "data") {
      variables$alpha_Mo = alpha
      Rs_Mo <- aggregate(data$Rs, as.yearmon(data$Date.daily, "%m/%y"), FUN = sum)
      R_W <- (1 - variables$alpha_Mo) * Rs_Mo
    }
 
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
      ET_type <- "Potential ET"
    } else if (est == "shallow lake") {
      ET_Mo.Monthly <- E_W
      ET_Mo.Average <- E_W.temp
      ET_type <- "Shallow Lake Evaporation"
    }

    ET.Daily <- NULL
    ET.Monthly <- ET_Mo.Monthly
    ET.Annual <- aggregate(ET.Monthly, floor(as.numeric(as.yearmon(data$Date.monthly, "%m/%y"))), FUN = sum)
    
    ET.MonthlyAve <- ET.AnnualAve <- NULL
    for (mon in min(as.POSIXlt(data$Date.monthly)$mon):max(as.POSIXlt(data$Date.monthly)$mon)){
      i = mon - min(as.POSIXlt(data$Date.monthly)$mon) + 1
      ET.MonthlyAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$mon== mon])
    }
    for (year in min(as.POSIXlt(data$Date.monthly)$year):max(as.POSIXlt(data$Date.monthly)$year)){
      i = year - min(as.POSIXlt(data$Date.monthly)$year) + 1
      ET.AnnualAve[i] <- mean(ET_Mo.Average[as.POSIXlt(data$Date.monthly)$year== year])
    }
    
    # Generate summary message for results
    ET_formulation <- "Morton CRWE"
    
    message(ET_formulation, " ", ET_type)
    message(variables$message1)
    message(variables$message6)

    results <- list(ET.Daily=ET.Daily, ET.Monthly=ET.Monthly, ET.Annual=ET.Annual, ET.MonthlyAve=ET.MonthlyAve, ET.AnnualAve=ET.AnnualAve, ET_formulation=ET_formulation, ET_type=ET_type, message1=variables$message1, message6=variables$message6)
    class(results) <- funname
    
    return(results)
  }
  # End of Morton's CRWE procedure
  
  
  #-------------------------------------------------------------------------------------

