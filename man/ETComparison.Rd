\name{ETComparison}
\alias{ETComparison}

\title{Plot esimtated evapotranspiration with climate variables
}
\description{
Produce comparison plots for results and statistics from different estimations produced by using different formulations and/or different input data. The number of different sets of results can be between 2 and 7. Plotting type can be selected among daily aggregation, monthly aggregation, annual aggregation, monthly average and annual average. For each type three comparison plots will be produced including time series, non-exceedance probability and box plot.
}
\usage{
ETComparison(results1, results2, results3 = NULL, results4 = NULL, results5 = NULL, 
results6 = NULL, results7 = NULL,labs, Sdate, Edate, type = "Monthly", ylim) 

}

\arguments{
  \item{results1}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor.
  }
  \item{results2}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor.
  }
  \item{results3}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor. The default is NULL if the user requires the comparison between only two sets of results.
  }
  \item{results4}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor. The default is NULL if the user requires the comparison among only three sets of results.
  }
  \item{results5}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor. The default is NULL if the user requires the comparison among only four sets of results.
  }
  \item{results6}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor. The default is NULL if the user requires the comparison among only five sets of results.
  }
  \item{results7}{
  A list named "results" which has been derived from function \code{ET..()} which can be of any class such as "Penman", "PenmanMonteith" and "PriestleyTaylor. The default is NULL if the user requires the comparison among only six sets of results.
  }
  \item{labs}{
  A character vector with the length equal to the number of sets of results to compare, defining the labels for the comparison plots}
  \item{Sdate}{
  Only used when type = 'Daily', 'Monthly' and 'Annual' to define the start date for the plotting windows, the default is the first day for the estimate evapotranspiration, but can be defined by user in the format "YYYY-MM-DD".
  }
  \item{Edate}{
  Only used when type = 'Daily', 'Monthly' and 'Annual' to define the end date for the plotting windows, the default is the last day for the estimate evapotranspiration, but can be defined by user in the format "YYYY-MM-DD".
  }
  \item{ylim}{
  A numeric vector of length 2 defining the lower and upper limit of the y-axis for plotting, the default is from 0 to 1.5 times of maximum value from the first set of result that is used to compare with others.
  }
  \item{type}{
  A character string indicating the type of plot produced, can be one of the following:
  "Daily" - comparison plots of estimated daily evapotranspiration;
  "Monthly" - comparison plots of monthly aggregated evapotranspiration;
  "Annual" - comparison plots of annually aggregated evapotranspiration;
  "MonthlyAve" - comparison plots of monthly averaged daily evapotranspiration;
  "AnnualAve" - comparison plots of annually averaged daily evapotranspiration.
  }
  }

\value{
Three plots are generated for each type of comparison plot selected, including:
1) time series plot of the estimated/aggregated/averaged values from each set of result;
2) non-exceedance plot of the distribution of estimated/aggregated/averaged values from each set of result; 
3) box plot of the distribution of estimated/aggregated/averaged values from each set of result.
  }

\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Penman" to call function
funname <- "Penman" 
class(data) <- funname 

# Call generic function ET(data, constants, ...) with class "Penman"
results_Penman <- ET(data, constants, solar="sunshine hours", wind=
"yes", windfunction_ver = "1948", alpha = 0.26, z0 = 0.02) 

# Set data class to be "PenmanMonteith" to call function
funname <- "PenmanMonteith" 
class(data) <- funname 

# Call generic function ET(data, constants, ...) with class "PenmanMonteith"
results_PenmanMonteith <- ET(data, constants, solar="sunshine hours", wind="yes", crop="short") 

# Plot the estimated Penman open-water evaporation against average temperature
ETComparison(results_Penman, results_PenmanMonteith, 
labs = c("Penman","PenmanMonteith"), Sdate = NULL, Edate = NULL, 
type = "Monthly", ylim=c(0,400))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{climate}
\keyword{forcing}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

