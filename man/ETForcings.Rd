\name{ETForcings}
\alias{ETForcings}

\title{Plot esimtated evapotranspiration with climate variables
}
\description{
Produce plot of daily, monthly and annual averaged estimated evapotranspiration with selected climate variables of the same time step.
}
\usage{
ETForcings(data, results, forcing)

}

\arguments{
  \item{data}{
  A list of data named \code{data} which must contain a component with the name of a climate variable that the estimated evapotranspiration should be plotted against, see \code{forcing}.
  }
  \item{results}{
  A list named \code{results} which has been derived from function \code{ET}.
  }
  \item{forcing}{
  A character string as the name of a climate variable that the estimated evapotranspiration should be plotted against, can be any of:\cr
  \emph{Tmax} - maximum temperature,\cr
  \emph{Tmin} - minimum temperature,\cr
  \emph{u2} - average wind speed at 2m,\cr
  \emph{uz} - average wind speed,\cr
  \emph{Rs} - solar radiation,\cr
  \emph{n} - daily sunshine hours,\cr
  \emph{Precip} - precipitation,\cr
  \emph{Epan} - Class-A pan evaporation,\cr
  \emph{RHmax} - maximum relative humidity,\cr
  \emph{RHmin} - minimum relative humidity,\cr
  \emph{Tdew} - average dew point temeprature.
  }
  }

\value{
Three plots are generated for the response of calculated evapotranspiration to each climate variable, including:\cr
1) daily evapotranspiration estimate vs. daily average temperature; \cr
2) monthly mean daily evaporationion estimate vs. monthly average temperature; \cr
3) annual mean daily evaporationion estimate vs. annual average temperature.
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
results <- ET.Penman(processeddata, constants, ts="daily", solar="sunshine hours",
wind="yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001, save.csv="no")

# Plot the estimated Penman open-water evaporation against average temperature
ETForcings(processeddata, results, forcing = "Tmax")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{climate}
\keyword{forcing}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

