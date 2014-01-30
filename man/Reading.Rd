\name{Reading}
\alias{Reading}

\title{Reading raw date and climate data
}
\description{
Load raw date and climate data, perform pre-processing, check for missing and error entries and then compile data list of daily time step.
}
\usage{
Reading(climatedata, constants, stopmissing)

}

\arguments{
  \item{climatedata}{
  A list of data named "climatedata" containing the raw data of date and climate variables. 
  The list must contain columns named as "Year", "Month" and "Day" to indicate the date.
  The climate variables will be of the following names while it is not compulsory to have all of them:
  "Tmax" - daily maximum temperature in degree Celcius,
  "Tmin" - daily minimum temperature in degree Celcius, 
  "Air.temperature.in.Degrees.C" - subdaily temperature in Degree Celcius,
  "Wind.speed.at.2m.measured.in.km.h" - subdaily wind speed at 2m from ground in km/h,
  "Wind.speed.measured.in.km.h" - subdaily wind speed measured at the height of wind instrument in km/h,
  "Solar.radiation.in.MJ/sqm/day" - daily incoming solar radiation,
  "n" - daily sunshine hours,
  "cloud.cover.in.oktas" - daily cloud cover in oktas,
  "precipitation.in.mm" - daily precipitation in mm,
  "Class.A.pan.evaporation.in.mm" - daily Class-A pan evaporation in mm,
  "RHmax" - daily maximum relative humidity in percentage,
  "RHmin" - daily minimum relative humidity in percentage,
  "Relative.humidity.in.percentage.." - subdaily relative humidity in percentage,
  "Dew.point.temperature.in.Degrees.C" - subdaily dew point temperature in degree Celcius,
  In order to determine which variables to include in "climatedata", please see \code{Evapotranspiration} for the specific data requirements for different formulations.
  }
  \item{constants}{
  A list named "constants" consists of constants required for data pre-processing which may contain the following items:
  "a_0", "b_0", "c_0", "d_0".
  These four constants which are constants required to calculate daily sunshine hours from daily cloud cover (Chiew and McMahon, 1991, described as a0, b0, c0, d0 in Equation S3.10 in McMahon et al., 2012), if such calculation will be performed these constants must be included in "constants".
  }
  \item{stopmissing}{
  A numeric value which represents the maximum percentage of missing data that the user can tolerate. It should be a number between 1 and 99. The percentage of missing data in each input variable is compared to this threshold, below which the missing data will be interpolated and otherwise the program will be terminated due to unsatisfactory data quality.
  }
  }

\value{
This function returns a list with all components of class \code{zoo} which have been processed from the raw data, including:
  \item{Date.daily}{
  A \code{zoo} object containing the date in daily step in the format of yyyy-mm-dd.
  }
  \item{Date.monthly}{
  A \code{zoo} object containing the date in daily step in the format of mmm-yyyy.
  }
  \item{J}{
  A \code{zoo} object containing the Julian Day for every day during the period that the data spans. 
  }
  \item{i}{
  A \code{zoo} object containing the month number for every day during the period that the data spans.
  }
  \item{ndays}{
  A \code{zoo} object containing the number of days for every month during the period that the data spans.
  }
  \item{Tmax}{
  A \code{zoo} object containing the daily maximum temperatures in degree Celcius.
  }
  \item{Tmin}{
  A \code{zoo} object containing the daily minimum temperatures in degree Celcius.
  }
  \item{u2}{
  A \code{zoo} object containing the daily wind speed at 2m from the ground in m/s.
  }
  \item{uz}{
  A \code{zoo} object containing the daily wind speed measured at the height of wind instrument in m/s.
  }
  \item{Rs}{
  A \code{zoo} object containing the daily solar radiation in MJ/m^2/day.
  }
  \item{n}{
  A \code{zoo} object containing the daily sunshine hours.
  }
  \item{Cd}{
  A \code{zoo} object containing the daily cloud cover in oktas.
  }
  \item{Precip}{
  A \code{zoo} object containing the daily precipitation in mm.
  }
  \item{Epan}{
  A \code{zoo} object containing the daily Class-A pan evaporation in mm.
  }
  \item{RHmax}{
  A \code{zoo} object containing the daily maximum relative humidity in percentage.
  }
  \item{RHmin}{
  A \code{zoo} object containing the daily minimum relative humidity in percentage.
  }
  \item{Tdew}{
  A \code{zoo} object containing the average daily dew temperatures in degree Celcius.
  }
  Note that the components might have value of \code{NULL} when the corresponding raw data cannot be found in "climatedata".
  }
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

CHIEW, F. H. & MCMAHON, T. A. 1991. THE APPLICABILITY OF MORTON'S AND PENMAN'S EVAPOTRANSPIRATION ESTIMATES IN RAINFALL-RUNOFF MODELING1. JAWRA Journal of the American Water Resources Association, 27, 611-620.
}
\author{
Danlu Guo
}
\seealso{
\code{evapotranspiration}
}
\examples{
# Reading climate data
data("climatedata")
data("constants")
data <- Reading(climatedata, constants, stopmissing = 40)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Reading}
\keyword{evapotranspiration}
\keyword{data pre-processing}
% __ONLY ONE__ keyword per line

