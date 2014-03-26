\name{climatedata}
\docType{data}
\alias{climatedata}
\title{Raw Climate Data Required for Calculating Evapotranspiration at the Kent Town station in Adelaide, Australia}
\description{
  This data set contains the raw climate data including the variables required for calculating evapotranspiration in function \code{ET} over the observation period between 1/3/2001 and 08/31/2004 at the Kent Town station in Adelaide, Australia.
}
\usage{climatedata}
\format{A data frame containing 10240 obserations of 13 objects:
"Station.Number" - weather station number,
"Year" - year of record,
"Month" - month of record,
"Day" - day of record,
"Hour" - hour of record,
"Julian" - Julian day of record,
"Temp.subdaily" - subdaily temperature data in degree Celcius,
"Tdew.subdaily" - subdaily dew point temperature data in degree Celcius,
"RH.subdaily" - subdaily relative humidity data in degree Celcius,
"n.daily" - daily sunshine hour data in hours,
"uz.subdaily" - subdaily wind speed data in kilometres per hour,
"Tmin.daily" - daily maximum temperature data in degree Celcius,
"Tmax.daily" - daily minimum temperature data in degree Celcius.
}
\source{Bureau of Meteorology, Kent Town, Adelaide, Australia}
\keyword{datasets}
