\name{ReadInputs}
\alias{ReadInputs}

\title{ReadInputs raw date and climate data
}
\description{
Load raw date and climate data, perform pre-processing, check for missing and error entries and then compile data list of daily time step.
}
\usage{
ReadInputs(varnames, climatedata, constants, stopmissing, timestep = "daily", 
                        interp_missing_days = FALSE, 
                        interp_missing_entries = FALSE, 
                        interp_abnormal = FALSE, 
                        missing_method = NULL, 
                        abnormal_method = NULL,
                        message = "yes")
}

\arguments{
  \item{varnames}{
  A character vector with length equals to the number of climate variables to be processed. Can include any element from:
  \emph{Tmax}, \emph{Tmin}, \emph{Temp}, \emph{Tdew},
  \emph{RHmax}, \emph{RHmin}, \emph{RH},
  \emph{Rs}, \emph{n}, \emph{Cd}, \emph{Precip},
  \emph{uz}, \emph{u2}, \emph{Epan}, \emph{va}, \emph{vs}.\cr
  
  Each variable is detailed as below:
  
  \emph{Tmax} - daily maximum temperature in degree Celcius,\cr
\emph{Tmin} - daily minimum temperature in degree Celcius,\cr
\emph{Temp} - subdaily temperature in degree Celcius\cr
\emph{Tdew} - dew point temperature in degree Celcius, either daily or subdaily accepted,\cr
\emph{RHmax} - daily maximum relative humidity in percentage,\cr
\emph{RHmin} - daily minimum relative humidity in percentage,\cr
\emph{RH} - subdaily relative humidity in degree Celcius,\cr
\emph{Rs} - incoming solar radiation in Megajoules per square metres per day, either daily or subdaily accepted,\cr
\emph{n} - daily sunshine hour in hours,\cr
\emph{Cd} - daily cloud cover in oktas,\cr
\emph{Precip} - precipitation in millimitres, either daily or subdaily accepted,\cr
\emph{u2} - wind speed measured at 2 metres from the ground surface in metres per second, either daily or subdaily accepted,\cr
\emph{uz} - wind speed in metres per second, either daily or subdaily accepted,\cr
\emph{Epan} - daily Class-A pan evaporation in millimitres,\cr
\emph{va} - average vapour pressure in KPa, either daily or subdaily accepted,\cr
\emph{vs} - saturated vapour pressure in KPa, either daily or subdaily accepted.\cr

  }
  \item{climatedata}{A data frame named "climatedata" containing the raw data of date and climate variables. \cr
  The data frame must have objects named as \emph{Year}, \emph{Month} and \emph{Day} to indicate the date.\cr
  The climate variables to include should be consistent with \emph{varnames}. \cr
  In order to determine which variables are needed for ET estimation, please see \code{\link{ET}} for the specific data requirements for different formulations.
  
  }
  \item{timestep}{
  Should be either \code{daily} or \code{subdaily} to specify the time step of raw climate data used.
  }
  \item{constants}{
  A list named "constants" consists of constants required for data pre-processing which may contain the following items:\cr
  \emph{a_0}, \emph{b_0}, \emph{c_0}, \emph{d_0}.\cr
  These four constants which are constants required to calculate daily sunshine hours from daily cloud cover (see Equation S3.10 in McMahon et al., 2012) - if the user requires such calculation these constants must be included in "constants".\cr
  The suggested values for various Australian locations are presented in Chiew and McMahon (1991), in which the four constants are named as a0, b0, c0, d0.
  }
  \item{stopmissing}{
  A numeric vector of length 3:\cr
  - the first value represents the maximum percentage of missing data that the user can tolerate;\cr
  - the second value represents the maximum percentage of the duration of missing data to the total data duaration that the user can toleratre;\cr
  - the third value represents the maximum percentage of missing days (within the date data, as a fraction of the total number of days) that the user can tolerate.\cr
  All values should be numbers between 1 and 99. \cr
  The percentages of the number and duration of missing data in the date data and each input variable are compared to the corresponding threshold; if any of the threshold is exceeded the program will be terminated due to unsatisfactory data quality.
  }
  \item{interp_missing_days}{
  \code{T} or \code{F}, indicating if missing days (within the date data) should be interpolated, with a default of \code{F} which assigns \code{NA} to data at the missing days.\cr
  }
  \item{interp_missing_entries}{
  \code{T} or \code{F}, indicating if missing data entries within individual climate variables should be interpolated, with a default of \code{F} which assigns \code{NA} to the missing entries.
  }
  \item{interp_abnormal}{
  \code{T} or \code{F}, if abnormal values within individual climate variables should be interpolated, with a default of \code{F} which leaves the abnormal values unchanged.\cr
  Abnormal values are defined differently according to the input variable, as following:\cr
  \cr
  - \emph{Tmax} > 100 or < -50 degree Celcius\cr
  - \emph{Tmin} > Tmax or < -50 degree Celcius \cr
  - \emph{Temp} > 100 or < -50 degree Celcius\cr
  - \emph{Tdew} > 100 or < -50 degree Celcius\cr
  - \emph{RHmax} > 100 or < 0 per cent\cr
  - \emph{RHmin}  > RHmax or < 0 per cent\cr
  - \emph{RH} > 100 per cent\cr 
  - \emph{Rs} < 0 MJ.m^2\cr
  - \emph{n} < 0 hour\cr
  - \emph{Cd}  < 0 Okta\cr
  - \emph{Precip}  < 0 mm\cr
  - \emph{uz}	< 0 m/s\cr
  - \emph{u2}	< 0 m/s\cr
  - \emph{Epan}	< 0 mm\cr
  - \emph{vs}	< 0 KPa\cr
  - \emph{va}	< 0 KPa\cr
  
  }
  \item{missing_method}{
  A character string for the name of the interpolated methods chosen for filling in missing days and missing data entries. Can be either: \cr
  \cr
  \code{monthly average} - replacement with same-month average (adapted from Narapusetty et al., 2009); \cr
  \code{seasonal average} - replacement with same-season average (adapted from Narapusetty et al., 2009);\cr
  \code{DoY average} - replacement with same day-of-the-year average (Narapusetty et al., 2009);\cr
  \code{neighbouring average} - interpolation between the two bounding values, which is only suitable for time increments in which values are available at adjacent increments (McMahon et al., 2013). When there is more than one consecutive missing entry, this interpolation fails, with a warning given.
  
  }
  \item{abnormal_method }{
  A character string for the name of the interpolated methods chosen for replacing data entries with abnormal values. Can be either: \cr
  \cr
  \code{monthly average} - replacement with same-month average (adapted from Narapusetty et al., 2009); \cr
  \code{seasonal average} - replacement with same-season average (adapted from Narapusetty et al., 2009);\cr
  \code{DoY average} - replacement with same day-of-the-year average (Narapusetty et al., 2009);\cr
  \code{neighbouring average} - interpolation between the two bounding values, which is only suitable for time increments in which non-abnormal values are available at adjacent increments (McMahon et al., 2013). When there is more than one consecutive abnormal entry, this interpolation fails, with a warning given.
  
  }
  \item{message }{"yes" or "no" indicating whether checking messages should be printed on screen.
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
  Note that the components might have value of \code{NULL} when the corresponding input variable cannot be found in the raw data (i.e. "climatedata").
  }
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Chiew, F. H. & McMahon, T. A. 1991. \emph{The applicability of Morton's and Penman's evapotranspiration estimates in rainfall-runoff modeling1}. JAWRA Journal of the American Water Resources Association, 27, 611-620.

Narapusetty, B., DelSole, T.Tippett, M.K. 2009, \emph{Optimal Estimation of the Climatological Mean}. Journal of Climate, vol. 22, no. 18, pp. 4845-4859.

}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{climatedata}},\code{\link{data}}
}
\examples{
# ReadInputs climate data
data("climatedata")
data("constants")
data <- ReadInputs(varnames = c("Temp","Tdew","n","RH","uz"),
                   climatedata, 
                   constants, 
                   stopmissing=c(10,10,3),
                   timestep = "subdaily",
                   interp_missing_days = FALSE, 
                   interp_missing_entries = FALSE, 
                   interp_abnormal = FALSE, 
                   missing_method = NULL, 
                   abnormal_method = NULL)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ReadInputs}
\keyword{evapotranspiration}
\keyword{data pre-processing}
% __ONLY ONE__ keyword per line
