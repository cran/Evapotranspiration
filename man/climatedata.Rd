\name{climatedata}
\docType{data}
\alias{climatedata}
\title{Raw Climate Data Required for Calculating Evapotranspiration}
\description{
  This data set contains the raw climate data including the variables required for calculating evapotranspiration in function \code{\link{ET}} over the observation period between 1/3/2001 and 08/31/2004 at the Kent Town station in Adelaide, Australia.
}
\usage{climatedata}
\format{A data frame containing 10240 obserations of 13 objects:\cr
\cr
\emph{Station.Number} - weather station number,\cr
\emph{Year} - year of record,\cr
\emph{Month} - month of record,\cr
\emph{Day} - day of record,\cr
\emph{Hour} - hour of record,\cr
\emph{Julian} - Julian day of record,\cr
\emph{Temp.subdaily} - subdaily temperature data in degree Celcius,\cr
\emph{Tdew.subdaily} - subdaily dew point temperature data in degree Celcius,\cr
\emph{RH.subdaily} - subdaily relative humidity data in degree Celcius,\cr
\emph{n.daily} - daily sunshine hour data in hours,\cr
\emph{uz.subdaily} - subdaily wind speed data in meter per second,\cr
\emph{Tmin.daily} - daily maximum temperature data in degree Celcius,\cr
\emph{Tmax.daily} - daily minimum temperature data in degree Celcius.
}
\source{Bureau of Meteorology, Kent Town, Adelaide, Australia}
\keyword{datasets}
