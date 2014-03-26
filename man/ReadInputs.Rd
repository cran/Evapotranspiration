\name{ReadInputs}
\alias{ReadInputs}

\title{ReadInputs raw date and climate data
}
\description{
Load raw date and climate data, perform pre-processing, check for missing and error entries and then compile data list of daily time step.
}
\usage{
ReadInputs(climatedata, constants, stopmissing)

}

\arguments{
  \item{climatedata}{
  A data frame named "climatedata" containing the raw data of date and climate variables. 
  The data frame contain objects named as "Year", "Month" and "Day" to indicate the date.
  The climate variables will be of the following names while it is not compulsory to have all of them:
  
  "Tmax.daily" - daily maximum temperature in degree Celcius,
  "Tmin.daily" - daily minimum temperature in degree Celcius, 
  "Temp.subdaily" - subdaily temperature in degree Celcius,
  "Tdew.subdaily" - subdaily dew point temperature in degree Celcius,
  "RHmax.daily" - daily maximum relative humidity in percentage,
  "RHmin.daily" - daily minimum relative humidity in percentage,
  "RH.subdaily" - subdaily relative humidity in degree Celcius,
  "Rs.daily" - daily incoming solar radiation in Millijoules per square metres per day,
  "n.daily" - daily sunshine hour in hours,
  "Cd.daily" - daily cloud cover in oktas,
  "Precip.daily" - daily precipitation in millimitres,
  "u2.subdaily" - subdaily wind speed measured at 2 metres from the ground surface in kilometres per hour,
  "uz.subdaily" - subdaily wind speed in kilometres per hour,
  "Epan.daily" - daily Class-A pan evaporation in millimitres,
  "Vp.subdaily" - subdaily vapour pressure in hectopascal.
  
  In order to determine which variables to include in "climatedata", please see \code{Evapotranspiration} for the specific data requirements for different formulations.
  }
  \item{constants}{
  A list named "constants" consists of constants required for data pre-processing which may contain the following items:
  "a_0", "b_0", "c_0", "d_0".
  These four constants which are constants required to calculate daily sunshine hours from daily cloud cover (Chiew and McMahon, 1991, described as a0, b0, c0, d0 in Equation S3.10 in McMahon et al., 2012), if such calculation will be performed these constants must be included in "constants".
  }
  \item{stopmissing}{
  A numeric vector of length 2. The first value represents the maximum percentage of missing data that the user can tolerate; the second value represents the maximum percentage of the duration of missing data to the total data duaration that the user can toleratre. Both values should be numbers between 1 and 99. The percentages of the number and duration of missing data in each input variable are compared to the corresponding threshold, below which the missing data will be interpolated and otherwise the program will be terminated due to unsatisfactory data quality.
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
# ReadInputs climate data
data("climatedata")
data("constants")
data <- ReadInputs(climatedata, constants, stopmissing = c(10,3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ReadInputs}
\keyword{evapotranspiration}
\keyword{data pre-processing}
% __ONLY ONE__ keyword per line

