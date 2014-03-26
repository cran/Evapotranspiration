ETPlot <- function(results, type = "Aggregation", OBS, OBSplot = FALSE,  Sdate = time(results$PET.Daily)[1], Edate = time(results$PET.Daily)[length(results$PET.Daily)]) { # plot estimations and observations
  
  # Aggregation plot (default)
  if (type == "Aggregation") {
    if (is.null(Sdate)) {
      Sdate = time(results$PET.Daily)[1]
    }
    if (is.null(Edate)) {
      Edate = time(results$PET.Daily)[length(results$PET.Daily)]
    }
    # Plots - Aggregation
    par(ask=FALSE)
    plot.new()
    x.Date <- as.Date(time(results$PET.Daily[which(time(results$PET.Daily) == as.Date(Sdate)):which(time(results$PET.Daily) == as.Date(Edate))]))
    if (!is.null(results$PET.Daily)) {
      plot(results$PET.Daily[which(as.Date(time(results$PET.Daily)) == x.Date[1]) : which(as.Date(time(results$PET.Daily)) == x.Date[length(x.Date)])], main = paste("Daily", results$PET_formulation, results$PET_type), xlab = "Year", ylab = list(c(results$PET_type, "mm/day")))
      if (OBSplot == TRUE) {
        if (!is.null(OBS)) {
          if (!is.null(OBS$E_obs.Daily)) {
            lines(OBS$E_obs.Daily[which(as.Date(time(results$PET.Daily)) == x.Date[1]) : which(as.Date(time(results$PET.Daily)) == x.Date[length(x.Date)])], main = "Observed evaporation", type = "o", pch = ".", col = "RED", xlab = "Year", ylab = "Observed evaporation mm/day")
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
    
    plot(results$PET.Monthly[which(as.yearmon(time(results$PET.Monthly)) == as.yearmon(x.Date[1])):which(as.yearmon(time(results$PET.Monthly)) == as.yearmon(x.Date[length(x.Date)]))], ylim = c(0,max(results$PET.Monthly)*1.5), main = paste("Monthly", results$PET_formulation, results$PET_type), type = "o", pch = ".", xlab = "Year", ylab = list(c(results$PET_type, "mm/month")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.Monthly[which(as.yearmon(time(results$PET.Monthly)) == as.yearmon(x.Date[1])):which(as.yearmon(time(results$PET.Monthly)) == as.yearmon(x.Date[length(x.Date)]))], main = "Observed evaporation", type = "o", pch = ".", ylim = c(0,400), col = "RED", xlab = "Year", ylab = "Observed evaporation mm/month")
        legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    plot(results$PET.Annual[which(floor(as.numeric(as.yearmon(time(results$PET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[1])))):which(floor(as.numeric(as.yearmon(time(results$PET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[length(x.Date)]))))],  main = paste("Annual", results$PET_formulation, results$PET_type), type = "o", xlab = "Year", ylim = c(0,3000), ylab = list(c(results$PET_type, "mm/year")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.Annual[which(floor(as.numeric(as.yearmon(time(results$PET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[1])))):which(floor(as.numeric(as.yearmon(time(results$PET.Annual)))) == floor(as.numeric(as.yearmon(x.Date[length(x.Date)]))))], main = "Observed evaporation", type = "o", col = "RED", ylim = c(0,3000), xlab = "Year", ylab = "Observed evaporation mm/year")
        legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=FALSE)
    paste("Completed plotting aggregations for results calculated by", results$PET_formulation, "formulation")
  }
  
  # Aggregation plot (default)
  if (type == "Average") {
    Sdate = time(results$PET.Daily)[1]
    Edate = time(results$PET.Daily)[length(results$PET.Daily)]
    
    par(ask=FALSE)
    plot.new()
    # Plots - Average
    plot(results$PET.MonthlyAve~unique(1:12), main = paste("Monthly average", results$PET_formulation, results$PET_type), type = "o", ylim = c(0, max(results$PET.MonthlyAve)*1.5), xlab = "Month", ylab = list(c("Monthly average", results$PET_type, "mm/day")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.MonthlyAve~unique(1:12), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Month", ylab = "Monthly average observed evaporation mm")
        legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=TRUE)
    
    plot(results$PET.AnnualAve~unique(as.POSIXlt(data$Date.daily)$year + 1900), main = paste("Annual average", results$PET_formulation, results$PET_type), type = "o", ylim = c(0, max(results$PET.AnnualAve)*1.5), xlab = "Year", ylab = list(c("Annual average", results$PET_type, "mm/day")))
    if (OBSplot == TRUE) {
      if (!is.null(OBS)) {
        lines(OBS$E_obs.AnnualAve~unique(as.POSIXlt(OBS$Date.OBS)$year + 1900), main = "Observed evaporation", type = "o", ylim = c(0,10), col = "RED", xlab = "Year", ylab = "Annual average observed evaporation mm")
        legend("topright", inset = .05, c(paste(results$PET_formulation,results$PET_type),"Observed Class-A pan evaporation"), cex = 0.8, col = c("BLACK","RED"), lty = 1)
      } else {
        warning("No observed data available for plotting, only the estimated values are plotted")
      }
    }
    
    par(ask=FALSE)
    paste("Completed plotting averages for results calculated by", results$PET_formulation, "formulation")
  }
}

#-------------------------------------------------------------------------------------

ETComparison <- function(results1, results2, results3 = NULL, results4 = NULL, results5 = NULL, results6 = NULL, results7 = NULL, type = "Monthly", ymin = NULL, ymax = NULL) { # plot up to 7 estimations
  # check number of sets of results to compare and compile list of legend text
  
  if (exists("results1")==F | exists("results2")==F) {
    stop("Please provide at least results1 and results2 for producing comparison plot")
  }
  if (is.null(results1) | is.null(results2)) {
    stop("Please provide at least results1 and results2 for producing comparison plot")
  }
  
  Ncomp <- 2
  legtext <- vector(mode = "character", length = Ncomp)
  legtext[1] <- paste(results1$PET_formulation, results1$PET_type)
  legtext[2] <- paste(results2$PET_formulation, results2$PET_type)
  
  if (!is.null(results3)) {
    Ncomp <- 3
    length(legtext) <- 3
    legtext[3] <- paste(results3$PET_formulation, results3$PET_type)
  }
  if (!is.null(results4)) {
    Ncomp <- 4
    length(legtext) <- 4
    legtext[4] <- paste(results4$PET_formulation, results4$PET_type)
  }
  if (!is.null(results5)) {
    Ncomp <- 5
    length(legtext) <- 5
    legtext[5] <- paste(results5$PET_formulation, results5$PET_type)
  }
  if (!is.null(results6)) {
    Ncomp <- 6
    length(legtext) <- 6
    legtext[6] <- paste(results6$PET_formulation, results6$PET_type)
  }
  if (!is.null(results7)) {
    Ncomp <- 7
    length(legtext) <- 7
    legtext[7] <- paste(results7$PET_formulation, results7$PET_type)
  }
  
  # Time series
  if (type == "Daily") {
    # Daily
    par(ask=FALSE)
    if (is.null(results1$PET.Daily) | is.null(results2$PET.Daily)) {
      stop("Unable to compare results because estimations of daily ET are not available in both results1 and results2")
    }
    
    comp <- matrix(nrow=length(results1$PET.Daily),ncol=Ncomp+1)
    comp[,1] <- as.Date(time(results1$PET.Daily))
    comp[,2] <- results1$PET.Daily
    comp[,3] <- results2$PET.Daily
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$PET.Daily)) {
          comp[,4] <- results3$PET.Daily
        } else {
          warning("No estimated daily ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$PET.Daily)) {
          comp[,5] <- results4$PET.Daily
        } else {
          warning("No estimated daily ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$PET.Daily)) {
          comp[,6] <- results5$PET.Daily
        } else {
          warning("No estimated daily ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$PET.Daily)) {
          comp[,7] <- results6$PET.Daily
        } else {
          warning("No estimated daily ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$PET.Daily)) {
          comp[,8] <- results7$PET.Daily
        } else {
          warning("No estimated daily ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ymin)) {
      ymin <- 0
    }
    if (is.null(ymax)) {
      ymax <- max(results1$PET.Daily)*1.5
    }
    
    plot(comp[,2]~as.Date(comp[,1]), type = "l", main = "Daily estimations of ET",  ylim = c(ymin,ymax), col = 2, xlab = "Year", ylab = "Estimated ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~as.Date(comp[,1]), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$PET.Daily),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
        sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the daily ET", xlim = c(0,1), ylim = c(ymin,ymax), col = 2, xlab = "Non-exceedance probability", ylab = "Estimated daily ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the daily estimations of ET", xlab = "ET estmation method", ylab = "Distribution of daily PET in mm", show.names = FALSE, pars = list(ylim = c(ymin,ymax), cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
  }

  if (type == "Monthly") {
    # Monthly
    par(ask=FALSE)
    if (is.null(results1$PET.Monthly) | is.null(results2$PET.Monthly)) {
      stop("Unable to compare results because estimations of monthly ET are not available in both results1 and results2")
    }
    
    comp <- matrix(nrow=length(results1$PET.Monthly),ncol=Ncomp+1)
    comp[,1] <- as.yearmon(time(results1$PET.Monthly))
    comp[,2] <- results1$PET.Monthly
    comp[,3] <- results2$PET.Monthly
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$PET.Monthly)) {
          comp[,4] <- results3$PET.Monthly
        } else {
          warning("No estimated monthly ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$PET.Monthly)) {
          comp[,5] <- results4$PET.Monthly
        } else {
          warning("No estimated monthly ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$PET.Monthly)) {
          comp[,6] <- results5$PET.Monthly
        } else {
          warning("No estimated monthly ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$PET.Monthly)) {
          comp[,7] <- results6$PET.Monthly
        } else {
          warning("No estimated monthly ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$PET.Monthly)) {
          comp[,8] <- results7$PET.Monthly
        } else {
          warning("No estimated monthly ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ymin)) {
      ymin <- 0
    }
    if (is.null(ymax)) {
      ymax <- max(results1$PET.Monthly)*1.5
    }
    
    plot(comp[,2]~as.yearmon(comp[,1]), type = "l", main = "Monthly estimations of ET",  ylim = c(ymin,ymax), col = 2, xlab = "Year", ylab = "Monthly aggregations of estimated ET mm/month", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~as.yearmon(comp[,1]), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$PET.Monthly),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the monthly ET", xlim = c(0,1), ylim = c(ymin,ymax), col = 2, xlab = "Non-exceedance probability", ylab = "Monthly aggregations of estimated ET mm/month", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), main = "Box plot of the monthly estimations of ET", boxwex = 0.5, range = 0, xlab = "ET estmation method", ylab = "Distribution of monthly PET in mm", show.names = FALSE, pars=list(ylim = c(ymin,ymax), cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
  }
  
  if (type == "Annual") {
    # Annual
    par(ask=FALSE)
    if (is.null(results1$PET.Annual) | is.null(results2$PET.Annual)) {
      stop("Unable to compare results because estimations of annual ET are not available in both results1 and results2")
    }
    
    comp <- matrix(nrow=length(results1$PET.Annual),ncol=Ncomp+1)
    comp[,1] <- as.yearmon(time(results1$PET.Annual))
    comp[,2] <- results1$PET.Annual
    comp[,3] <- results2$PET.Annual
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$PET.Annual)) {
          comp[,4] <- results3$PET.Annual
        } else {
          warning("No estimated annual ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$PET.Annual)) {
          comp[,5] <- results4$PET.Annual
        } else {
          warning("No estimated annual ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$PET.Annual)) {
          comp[,6] <- results5$PET.Annual
        } else {
          warning("No estimated annual ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$PET.Annual)) {
          comp[,7] <- results6$PET.Annual
        } else {
          warning("No estimated annual ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$PET.Annual)) {
          comp[,8] <- results7$PET.Annual
        } else {
          warning("No estimated annual ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ymin)) {
      ymin <- 0
    }
    if (is.null(ymax)) {
      ymax <- max(results1$PET.Annual)*1.5
    }
    
    plot(comp[,2]~floor(as.numeric(as.yearmon(comp[,1]))), type = "l", main = "Annual estimations of ET",  ylim = c(ymin,ymax), col = 2, xlab = "Year", ylab = "Annual aggregations of estimated ET mm/year", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~floor(as.numeric(as.yearmon(comp[,1]))), type = "l", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$PET.Annual),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the annual ET", xlim = c(0,1), ylim = c(ymin,ymax), col = 2, xlab = "Non-exceedance probability", ylab = "Annual aggregations of estimated ET mm/year", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the annual estimations of ET", xlab = "ET estmation method", ylab = "Distribution of annual PET in mm", show.names = FALSE, pars = list(ylim = c(ymin,ymax), cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
  }
  
  if (type == "MonthlyAve") {
    # Monthly average
    par(ask=FALSE)
    if (is.null(results1$PET.MonthlyAve) | is.null(results2$PET.MonthlyAve)) {
      stop("Unable to compare results because estimations of monthly averaged ET are not available in both results1 and results2")
    }
    
    comp <- matrix(nrow=length(results1$PET.MonthlyAve),ncol=Ncomp+1)
    comp[,1] <- as.yearmon(time(results1$PET.MonthlyAve))
    comp[,2] <- results1$PET.MonthlyAve
    comp[,3] <- results2$PET.MonthlyAve
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$PET.MonthlyAve)) {
          comp[,4] <- results3$PET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$PET.MonthlyAve)) {
          comp[,5] <- results4$PET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$PET.MonthlyAve)) {
          comp[,6] <- results5$PET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$PET.MonthlyAve)) {
          comp[,7] <- results6$PET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$PET.MonthlyAve)) {
          comp[,8] <- results7$PET.MonthlyAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ymin)) {
      ymin <- 0
    }
    if (is.null(ymax)) {
      ymax <- max(results1$PET.MonthlyAve)*1.5
    }
    
    plot(comp[,2]~unique(1:12), type = "l", main = "Monthly average ET",  ylim = c(ymin,ymax), col = 2, xlab = "Month", ylab = "Monthly averages of estimated daily ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~unique(1:12), type = "o", pch = ".", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$PET.MonthlyAve),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the monthly average ET", xlim = c(0,1), ylim = c(ymin,ymax), col = 2, xlab = "Non-exceedance probability", ylab = "Monthly averages of estimated daily ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the monthly average ET", xlab = "ET estmation method", ylab = "Distribution of monthly averaged PET in mm/day", show.names = FALSE, pars = list(ylim = c(ymin,ymax), cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
  }
  
  if (type == "AnnualAve") {
    # Annual average
    par(ask=FALSE)
    if (is.null(results1$PET.AnnualAve) | is.null(results2$PET.AnnualAve)) {
      stop("Unable to compare results because estimations of monthly averaged ET are not available in both results1 and results2")
    }
    
    comp <- matrix(nrow=length(results1$PET.AnnualAve),ncol=Ncomp+1)
    comp[,1] <- as.yearmon(time(results1$PET.AnnualAve))
    comp[,2] <- results1$PET.AnnualAve
    comp[,3] <- results2$PET.AnnualAve
    if (Ncomp >= 3) {
      if (!is.null(results3)) {
        if (!is.null(results3$PET.AnnualAve)) {
          comp[,4] <- results3$PET.AnnualAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results4")
          comp[,4] <- NULL
        }       
      }
    }
    if (Ncomp >= 4) {
      if (!is.null(results4)) {
        if (!is.null(results4$PET.AnnualAve)) {
          comp[,5] <- results4$PET.AnnualAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results5")
          comp[,5] <- NULL
        }       
      }
    }
    if (Ncomp >= 5) {
      if (!is.null(results5)) {
        if (!is.null(results5$PET.AnnualAve)) {
          comp[,6] <- results5$PET.AnnualAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results6")
          comp[,6] <- NULL
        }       
      }
    }
    if (Ncomp >= 6) {
      if (!is.null(results6)) {
        if (!is.null(results6$PET.AnnualAve)) {
          comp[,7] <- results6$PET.AnnualAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results7")
          comp[,7] <- NULL
        }       
      }
    }
    if (Ncomp >= 7) {
      if (!is.null(results7)) {
        if (!is.null(results7$PET.AnnualAve)) {
          comp[,8] <- results7$PET.AnnualAve
        } else {
          warning("No estimated monthly averaged ET available for plotting from results8")
          comp[,8] <- NULL
        }       
      }
    }
    
    if (is.null(ymin)) {
      ymin <- 0
    }
    if (is.null(ymax)) {
      ymax <- max(results1$PET.AnnualAve)*1.5
    }
    
    plot(comp[,2]~floor(as.numeric(as.yearmon(time(results1$PET.Annual)))), type = "l", main = "Annually average ET",  ylim = c(ymin,ymax), col = 2, xlab = "Month", ylab = "Annual averages of estimated daily ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (i in 1:(Ncomp-1)) {
      lines(comp[,i+2]~floor(as.numeric(as.yearmon(time(results1$PET.Annual)))), type = "l", col = i+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Non-exceedance probability
    par(ask=TRUE)
    sortcomp <- matrix(nrow=length(results1$PET.AnnualAve),ncol=Ncomp+2)
    
    sortcomp[,1] <- comp[,1]
    for (j in 1:Ncomp) {
      sortcomp[,j+1] <- sort(comp[,j+1])
    }
    sortcomp[,Ncomp+2] <- seq(from = 1/nrow(sortcomp), to = 1, by = 1/nrow(sortcomp))
    
    plot(sortcomp[,2]~sortcomp[,Ncomp+2], type = "l", main = "Non-exceedance probability of the annual average ET", xlim = c(0,1), ylim = c(ymin,ymax), col = 2, xlab = "Non-exceedance probability", ylab = "Annual averages of estimated daily ET mm/day", cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1)))
    for (k in 1:(Ncomp-1)) {
      lines(sortcomp[,k+2]~sortcomp[,Ncomp+2], type = "o", pch = ".", col = k+2)
    }
    
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
    
    # Box plot
    par(ask=TRUE)
    colnames(comp)[2:(Ncomp+1)] <- legtext
    boxplot(comp[,2:(Ncomp+1)], col = c(seq(from=2,to=Ncomp+1)), boxwex = 0.5, range = 0, main = "Box plot of the annual average ET", xlab = "ET estmation method", ylab = "Distribution ofannually averaged PET in mm/day", show.names = FALSE, pars = list(ylim = c(ymin,ymax), cex.lab = 1.9, cex.main = 2.1, par(mar = c(5.1,4.7,4.7,2.1))))
    legend("topright", inset = .01, c(paste(legtext[1:Ncomp])), cex = (5.5/Ncomp), col = c(seq(from=2,to=Ncomp+1)), lty = 1)
  }
}


#-------------------------------------------------------------------------------------

ETForcings <- function(data, results, forcing) {
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