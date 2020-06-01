\name{climatedata}
\alias{climatedata}
\docType{data}
\title{Raw Climate Data Required for Calculating Evapotranspiration}
\description{
  This data set contains the raw climate data including the variables required for calculating evapotranspiration in function \code{\link{ET}} over the observation period between 1/3/2001 and 08/31/2004 at the Kent Town station in Adelaide, Australia.
}
\usage{data(climatedata)}
\format{A data frame containing 10240 obserations of 9 objects:\cr
\cr
\emph{Station.Number} - weather station number,\cr
\emph{Year} - year of record,\cr
\emph{Month} - month of record,\cr
\emph{Day} - day of record,\cr
\emph{Hour} - hour of record,\cr
\emph{Tdew} - subdaily dew point temperature in degree Celcius,\cr
\emph{RH} - subdaily relative humidity in percentage,\cr
\emph{Rs} - subdaily solar radiation in Megajoule per square meter,\cr
\emph{uz} - subdaily wind speed in meter per second.

}
\source{Bureau of Meteorology, Kent Town, Adelaide, Australia}
\keyword{datasets}
