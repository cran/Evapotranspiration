\name{ETComparison}
\alias{ETComparison}

\title{Compare esimtated evapotranspiration among multiple sets of result
}
\description{
Produce comparison plots for results and statistics from different estimations produced by using different formulations and/or different input data. The number of different sets of results can be between 2 and 7. Plotting type can be selected among daily aggregation, monthly aggregation, annual aggregation, monthly average and annual average. For each type three comparison plots will be produced including time series, non-exceedance probability and box plot.
}
\usage{
ETComparison(results1, results2, results3 = NULL, results4 = NULL, results5 = NULL,
results6 = NULL, results7 = NULL, labs, Sdate = NULL, Edate = NULL,
type = "Monthly", ylim = rep(NA,2))

}

\arguments{
  \item{results1}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}.
  }
  \item{results2}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}.
  }
  \item{results3}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}. \cr
  The default is NULL if the user requires the comparison between only two sets of results.
  }
  \item{results4}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}. \cr
  The default is NULL if the user requires the comparison among only three sets of results.
  }
  \item{results5}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}. \cr
  The default is NULL if the user requires the comparison among only four sets of results.
  }
  \item{results6}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}.\cr
  The default is NULL if the user requires the comparison among only five sets of results.
  }
  \item{results7}{
  A list named \code{results} which has been derived from function \code{ET} which can be from any model such as \code{Penman}, \code{Penman-Monteith} or \code{Priestley-Taylor}. \cr
  The default is NULL if the user requires the comparison among only six sets of results.
  }
  \item{labs}{
  A character vector with the length equal to the number of sets of results to compare, defining the labels for the comparison plots}
  \item{Sdate}{
  Only used when argument \code{type} is \code{Daily}, \code{Monthly} or \code{Annual} to define the start date for the plotting windows, which can be defined by user in the format \emph{YYYY-MM-DD}; if missing the default is the first day of data is used.
  }
  \item{Edate}{
  Only used when argument \code{type} is \code{Daily}, \code{Monthly} or \code{Annual} to define the end date for the plotting windows, which can be defined by user in the format \emph{YYYY-MM-DD}; if missing the default is the last day of data is used.
  }
  \item{ylim}{
  A numeric vector of length 2 defining the lower and upper limit of the y-axis for plotting, if missing the default is from 0 to 1.5 times of maximum value from the first set of result that is used to compare with others.
  }
  \item{type}{
  A character string indicating the type of plot produced, can be one of the following:\cr
  \code{Daily} - comparison plots of estimated daily evapotranspiration;\cr
  \code{Monthly} - comparison plots of monthly aggregated evapotranspiration;\cr
  \code{Annual} - comparison plots of annually aggregated evapotranspiration;\cr
  \code{MonthlyAve} - comparison plots of monthly averaged daily evapotranspiration;\cr
  \code{AnnualAve} - comparison plots of annually averaged daily evapotranspiration.
  }
  }

\value{
Three plots are generated for each type of comparison plot selected, including:\cr
1) time series plot of the estimated/aggregated/averaged values from each set of result;\cr
2) non-exceedance plot of the distribution of estimated/aggregated/averaged values from each set of result; \cr
3) box plot of the distribution of estimated/aggregated/averaged values from each set of result.\cr
  }

\author{
Danlu Guo
}
\seealso{
\code{\link{ETPlot}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.Penman under the generic function ET
results_Penman <- ET.Penman(processeddata, constants, ts="daily", solar="sunshine hours",
wind="yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001, save.csv="no")

# Call ET.PenmanMonteith under the generic function ET
results_PenmanMonteith <- ET.PenmanMonteith(processeddata, constants, ts="daily",
solar="sunshine hours", wind="yes", crop = "short", save.csv="no")

# Plot the estimated Penman open-water evaporation against average temperature
ETComparison(results_Penman, results_PenmanMonteith, type = "Monthly", ylim=c(0,400),
labs=c("Penman","PenmanMonteith"))
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{climate}
\keyword{forcing}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

