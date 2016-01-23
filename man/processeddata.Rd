\name{processeddata}
\docType{data}
\alias{processeddata}
\title{Processed Climate Data Required for Calculating Evapotranspiration}
\description{
  This data set contains the processed climate data including the variables required for calculating evapotranspiration in function \code{\link{ET}} over the observation period between 1/3/2001 and 31/8/2004 at the Kent Town station in Adelaide, Australia.
}
\usage{processeddata}
\format{A list containing 17 non-empty climate variables:\cr
\cr
\emph{Date.daily} - date in daily time step,\cr
\emph{Date.monthly} - date in monthly time step,\cr
\emph{J} - julian days,\cr
\emph{i} - month,\cr
\emph{Ndays} - days in month,\cr
\emph{Tmax} - daily maximum temperature in degree Celcius,\cr
\emph{Tmin} - daily minimum temperature in degree Celcius,\cr
\emph{uz} - daily wind speed in meters per second,\cr
\emph{n} - daily sunshine hour in hours,\cr
\emph{RHmax} - daily maximum relative humidity in percentage,\cr
\emph{RHmin} - daily minimum relative humidity in percentage,\cr
\emph{Tdew} - daily dew point temperature data in degree Celcius.
}
\source{Bureau of Meteorology, Kent Town, Adelaide, Australia}
\keyword{datasets}
