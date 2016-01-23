\name{ETPlot}
\alias{ETPlot}

\title{Plot the daily, monthly and annual aggregations of esimtated evapotranspiration
}
\description{
Produce plot of aggregated estimations of evapotranspiration in daily, monthly and annual steps, or averaged daily estimations in monthly or annual steps. 
}
\usage{
ETPlot(results, type = "Aggregation", OBS, OBSplot, Sdate = time(results$ET.Daily)[1], 
Edate = time(results$ET.Daily)[length(results$ET.Daily)])

}

\arguments{
  \item{results}{
  A list named \code{results} which has been derived from function \code{ET..()}.
  }
  \item{type}{
  A character string of either \code{Aggregation} or \code{Average} to indicate the type of plot required. The default is \code{Aggregation}. For aggregation plot the user can define the start and end date of plotting or by default using the calculation period for plotting. For average plot the plotting period equals to the calculation period.
  }
  \item{OBS}{
  A list named \code{OBS} which has been derived from function \code{ReadOBSEvaporation}.
  }
  \item{OBSplot}{
  Must be eith TRUE or FALSE. TRUE indicates that the observed evaporation will be plotted together with the estimations and FALSE indicates that the observations will not be shown on the plots.
  }
  \item{Sdate}{
  Only used when type = \code{Aggregation} to define the start date for the plotting windows, the default is the first day for the estimate evapotranspiration, but can be defined by user in the format \emph{YYYY-MM-DD}.
  }
  \item{Edate}{
  Only used when type = \code{Aggregation} to define the end date for the plotting windows, the default is the last day for the estimate evapotranspiration, but can be defined by user in the format \emph{YYYY-MM-DD}.
  }
  }

\value{
If argument \code{type} is \code{Aggregation}, three plots are displayed in the following order (the next one appears after pressing enter): \cr
1)  Daily evapotranspiration estimates; \cr
2)	Monthly evapotranspiration estimates aggregated from daily estimates; \cr
3)	Annual evapotranspiration estimates aggregated from daily estimates. \cr
 \cr
If argument \code{type} is \code{Average}, two plots are displayed in the following order \cr
1) Monthly averaged daily estimations of evapotranspiration; \cr
2) Annually averaged daily estiamtions of evapotranspriation. \cr
  }

\author{
Danlu Guo
}
\seealso{
\code{\link{ETComparison}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.Penman under the generic function ET
results <- ET.Penman(data, constants, ts="daily", solar="sunshine hours", 
wind="yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001) 


# Read evaporation data
data("E_OBS")
OBS <- ReadOBSEvaporation(E_OBS, data)

# Plot the aggregation of estimated Penman open-water evaporation with observed evaporation
ETPlot(results, type = "Aggregation", OBS, OBSplot = TRUE, Sdate = "2001-05-01", 
Edate = "2004-05-01")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

