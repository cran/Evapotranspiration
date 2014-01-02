\name{PlotEvapotranspiration}
\alias{PlotEvapotranspiration}

\title{Plot esimtated evapotranspiration
}
\description{
Produce plot of aggregated and averaged estimated evapotranspiration in daily, monthly and annual steps.
}
\usage{
PlotEvapotranspiration(results, OBS, OBSplot)

}

\arguments{
  \item{results}{
  A list named "results" which has been derived from function \code{Evapotranspiration..()}.
  }
  \item{OBS}{
  A list named "OBS" which has been derived from function \code{ReadOBSEvaporation()}.
  }
  \item{OBSplot}{
  Must be eith TRUE or FALSE. TRUE indicates that the observed evaporation will be plotted together with the estimations and FALSE indicates that the observations will not be shown on the plots.
  }
  }

\value{
A series of plots are displayed one by one in the following order (the next one appears after pressing enter):
1)  Daily evapotranspiration estimates;
2)	Monthly evapotranspiration estimates aggregated from daily estimates;
3)	Annual evapotranspiration estimates aggregated from daily estimates;
4)	Monthly mean daily evapotranspiration averaged from daily estimates;
5)	Annual mean daily evapotranspiration averaged from daily estimates.

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
results <- Evapotranspiration(data, constants, sunshine.hours="data", wind=
"yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001) 

# Read evaporation data
data("E_OBS")
OBS <- ReadOBSEvaporation(E_OBS, data)

# Plot the estimated Penman open-water evaporation with observed evaporation
PlotEvapotranspiration(results, OBS, OBSplot = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{plot}
\keyword{evapotranspiration}
% __ONLY ONE__ keyword per line

