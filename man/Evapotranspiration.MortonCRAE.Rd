\name{Evapotranspiration.MortonCRAE}
\alias{Evapotranspiration.MortonCRAE}

\title{Morton CRAE Formulation
}
\description{
Implementing the Morton CRAE formulation for estimating potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
}
\usage{
\method{Evapotranspiration}{MortonCRAE}(data, constants, est, solar, Tdew, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "MortonCRAE" which contains the following items (climate variables) required by Morton CRAE formulation:
  "Tmax", "Tmin", "Tdew", "n" or "Cd".
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Morton CRAE formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lat_rad" - latitude in radians,
  "PA" - annual precipitation in mm, required when precipitation data is not available,
  "epsilonMo"" - surface emissivity = 0.92 for CRAE model (Morton, 1986)
  "fz" - A constant in Morton's procedure = 28.0 Wm^-2.mbar^-1 for T >= 0 degree Celcius, and = 28.0*1.15 Wm^-2.mbar^-1 for T >= 0 degree Celcius for CRAE model (Morton, 1983)
  "b0" - a constants in Morton's procedure, = 1 for CRAE model (Morton, 1983)
  "b1" - a constant in Morton's procedure, = 14 for CRAE model (Morton, 1983)
  "b2 - a constant in Morton's procedure, = 1.2 for CRAE model (Morton, 1983)
  "gammaps" - Produce of Psychrometric constant and atmospheric pressure as sea level, = 0.66 mbar. degree Celcius^-1 for CRAE model for T >= 0 degree Celcius, = 0.66/1.15 mbar. degree Celcius^-1 for CRAE model for T < 0 degree Celcius (Morton, 1983),
  "lat" - latitude in degrees,
  "alphaMo" - a constant in Morton's procedure, = 17.27 when T >= 0 degree Celcius, = 21.88 when T < 0 degree Celcius (Morton, 1983),
  "betaMo" - a constant in Morton's procedure, = 237.3 degree Celcius  when T >= 0 degree Celcius,, = 2655. degree Celcius, when T < 0 degree Celcius (Morton, 1983),
  "sigmaMo" - Stefan-Boltzmann constant in Morton's procedure, = 5.67e-08 W.m^-2.K^-4 (Morton, 1983),
  "lambdaMo" - Latent heat of vaporisation in Morton's procedure, = 28.5W.day.kg^-1 when T >= 0 degree Celcius, = 28.5*1.15W.day.kg^-1 when T < 0 degree Celcius.
  }
  \item{solar}{
  Must be either "sunshine hours" or "cloud". "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data.
  }
  \item{est}{
  Must be either "potential ET", "wet areal ET" or "actual areal ET". "potential ET" proceeds to estimating potential evapotranspiration; "wet areal ET" proceeds to estimating wet-environmental areal evapotranspiration; "actual areal ET" proceeds to estimating actual areal evapotranspiraion
  }
  \item{Tdew}{
  Logical, must be either TRUE or FALSE, indicating if real data of dew point temperature is used for calculating the radiation in Morton's formulations, if TRUE the data will be used and if FALSE the dew point temperature will be calculated from data of daily vapour pressure.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The type of evapotranspiration calculated can be selected through argument \code{est}, please see "Arguments" for details.
  The alternative calculation options can be selected through argument \code{solar} and \code{Tdew}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{MortonCRAE}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Morton CRAE potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Morton CRAE potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Morton CRAE potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Morton CRAE potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Morton CRAE potential evapotranspiration, wet-environment areal evapotranspiration and actual areal evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "MortonCRAE".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is either "Potential Evapotranspiration", "Wet-environment Areal Evapotranspiration" and "Actual Areal Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message6}{
  A message to inform the users about if actual dew point temperature has been used in the calculations or alternative calculations has been performed without dew point temperature data.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

MORTON, F. I. 1983. Operational estimates of areal evapotranspiration and their significance to the science and practice of hydrology. Journal of Hydrology, 66, 1-76.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "MortonCRAE" to call function
funname <- "MortonCRAE" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "MortonCRAE"
results <- Evapotranspiration(data, constants, est = "potential ET", solar = 
"sunshine hours", Tdew = TRUE) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MortonCRAE}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
\keyword{wet-environment areal evapotranspiration}
\keyword{actual areal evapotranspiration}
% __ONLY ONE__ keyword per line

