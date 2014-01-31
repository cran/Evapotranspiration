\name{EvapotranspirationForcings}
\alias{EvapotranspirationForcings}

\title{Plot esimtated evapotranspiration with climate variables
}
\description{
Produce plot of daily, monthly and annual averaged estimated evapotranspiration with selected climate variables of the same time step.
}
\usage{
EvapotranspirationForcings(data, results, forcing)

}

\arguments{
  \item{data}{
  A list of data named "data" which must contain a component with the name of a climate variable that the estimated evapotranspiration should be plotted against, see \code{forcing}.
  }
  \item{results}{
  A list named "results" which has been derived from function \code{Evapotranspiration..()}.
  }
  \item{forcing}{
  A character string as the name of a climate variable that the estimated evapotranspiration should be plotted against, can be any of:
  "Tmax" - maximum temperature,
  "Tmin" - minimum temperature,
  "u2" - average wind speed at 2m,
  "uz" - average wind speed,
  "Rs" - solar radiation,
  "n" - daily sunshine hours,
  "Precip" - precipitation,
  "Epan" - Class-A pan evaporation,
  "RHmax" - maximum relative humidity,
  "RHmin" - minimum relative humidity,
  "Tdew" - average dew point temeprature.
  }
  }

\value{
Three plots are generated for the response of calculated evapotranspiration to each climate variable, including:
1) daily evapotranspiration estimate vs. daily average temperature; 
2) monthly mean daily evaporationion estimate vs. monthly average temperature; 
3) annual mean daily evaporationion estimate vs. annual average temperature.
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

# Call generic function Evapotranspiration(data, constants, ...) with class "Penman"
results <- Evapotranspiration(data, constants, solar="sunshine hours", wind=
"yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001) 
# Plot the estimated Penman open-water evaporation against average temperature
EvapotranspirationForcings(data, results, forcing = "Tmax")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{climate}
\keyword{forcing}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

