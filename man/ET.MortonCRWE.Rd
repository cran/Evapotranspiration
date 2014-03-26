\name{ET.MortonCRWE}
\alias{ET.MortonCRWE}

\title{Morton CRWE Formulation
}
\description{
Implementing the Morton CRWE formulation for estimating potential evapotranspiration or shallow lake evaporation.
}
\usage{
\method{ET}{MortonCRWE}(data, constants, est, solar, Tdew, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "MortonCRWE" which contains the following items (climate variables) required by Morton CRWE formulation:
  "Tmax", "Tmin", "Tdew", "n" or "Cd".
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Morton CRWE formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lat_rad" - latitude in radians,
  "PA" - annual precipitation in mm, required when precipitation data is not available,
  "epsilonMo"" - surface emissivity = 0.92 for CRWE model (Morton, 1983)
  "fz" - A constant in Morton's procedure = 25.0 Wm^-2.mbar^-1 for T >= 0 degree Celcius, and = 28.75 Wm^-2.mbar^-1 for T >= 0 degree Celcius for CRWE model (Morton, 1983)
  "b0" - a constants in Morton's procedure, = 1.12 for CRWE model (Morton, 1983)
  "b1" - a constant in Morton's procedure, = 13 for CRWE model (Morton, 1983)
  "b2 - a constant in Morton's procedure, = 1.12 for CRWE model (Morton, 1983)
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
  Must be either "potential" or "shallow lake". "potential" proceeds to estimating potential evapotranspiration; "shallow lake" proceeds to estimating shallow lake evaporation.
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
This function returns a object of class \code{MortonCRWE}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "MortonCRWE".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is either "Potential Evapotranspiration" or "Shallow Lake Evaporation".
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

MORTON, F. I. 1986. Practical Estimates of Lake Evaporation. Journal of Applied Meteorology, 25, 371-388.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "MortonCRWE" to call function
funname <- "MortonCRWE" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "MortonCRWE"
results <- ET(data, constants, est = "potential", solar = 
"sunshine hours", Tdew = TRUE) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MortonCRWE}
\keyword{shallow lake}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

