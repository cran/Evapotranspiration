\name{ReadOBSEvaporations}
\alias{ReadOBSEvaporation}

\title{Read Raw Data of Observed Evaporation from file
}
\description{
Load raw date and evaporation data and then compile data list of daily time step.
}
\usage{
ReadOBSEvaporation(E_OBS, data)

}

\arguments{
  \item{E_OBS}{
  A list of evaporation data named \code{E_OBS} which must contain the following columns:\cr
  \code{Year}, \code{Month}, \code{Day} as the date and,\cr
  \code{EVAP.Obs} as the observed evaporation in mm.\cr
  The observations can be of daily and monthly time steps and must match with the corresponding dates recorded.
  }
  \item{data}{
  A list of data named \code{data} which contains data of climate variables over the same period as the evaporation data
  }
  }

\value{
This function returns a list with all components of class \code{zoo} which have been processed from the raw data, including:
  \item{Date.OBS}{
  A \code{zoo} object containing the date data with time step consistent with the raw evaporation data in \code{E_OBS}.
  }
  \item{E_obs.Daily}{
  A \code{zoo} object containing the daily evaporation data.
  }
  \item{E_obs.Monthly}{
  A \code{zoo} object containing the monthly aggregaated observed evaporation in mm.
  }
  \item{E_obs.Annual}{
  A \code{zoo} object containing the annually aggregated observed evaporation in mm.
  }
  \item{E_obs.MonthlyAve}{
  A \code{zoo} object containing the monthly averaged daily evaporation from observation in mm/day.
  }
  \item{E_obs.AnnualAve}{
  A \code{zoo} object containing the annually average daily evaporation from observation in mm/day.
  }
  Note that the components might have value of \code{NULL} when the corresponding raw data cannot be found in \code{E_OBS}.
  }

\author{
Danlu Guo
}

\examples{
# Get the time period from "data"
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Reading obsevations of evaporation within specified time period
data("E_OBS")
OBS <- ReadOBSEvaporation(E_OBS, processeddata)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reading}
\keyword{evapotranspiration}
\keyword{observation}
% __ONLY ONE__ keyword per line

