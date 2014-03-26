\name{ET.PenmanMonteith}
\alias{ET.PenmanMonteith}

\title{Penman-Monteith Formulation
}
\description{
Implementing the Penman-Monteith formulation (including the method for FAO-56 hypothetical short grass and the method for ASCE-EWRI Standardised crop) for estimating reference crop evapotranspiration
}
\usage{
\method{ET}{PenmanMonteith}(data, constants, solar, wind, crop, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "PenmanMonteith" which contains the following items (climate variables) required by Penman-Monteith formulation:
  "Tmax", "Tmin", "RHmax", "RHmin", "Rs, "n" or "Cd", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penman-Monteith formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "z" - height of wind instrument in m,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,
  "G" - soil heat flux in MJ.m^-2.day^-1, = 0 when using daily time step.
  }
  \item{solar}{
  Must be either "data", "sunshine hours", "cloud" or "monthly precipitation". "data" indicates that solar radiation data is to be used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{wind}{
  Must be either "yes" or "no". "yes" indicates that the calculation will use real data of wind speed; "no" indicates that the alternative calculation without using wind data will be used in Penman formulation (Valiantzas 2006, Equation33).
  }
  \item{crop}{
  Must be either "short" or "tall". "short" indicates that the method for FAO-56 hypothetical short grass will be applied (Allen et al., 1998, Equation 6); "tall" indicates that the method for  ASCE-EWRI Standardised crop will be applied (ASCE, 2005, Equation 1, Table 1)
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through arguments \code{solar} and \code{wind}, please see "Arguments" for details.
  User-defined evaporative surface is allowed through arguments \code{crop}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{PenmanMonteith}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Penman-Monteith rerference crop evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Penman-Monteith rerference crop evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Penman-Monteith rerference crop evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Penman-Monteith rerference crop evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Penman-Monteith rerference crop evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to either "Penman-Monteith FAO56" or "Penman-Monteith ASCE-EWRI Standardised".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message2}{
  A message to inform the users about if actual wind data has been used in the calculations or alternative calculations has been performed without wind data.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

ALLEN, R. G., PEREIRA, L. S., RAES, D. & SMITH, M. 1998. Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300, 6541.

ALLEN, R. G. 2005. The ASCE standardized reference evapotranspiration equation, Amer Society of Civil Engineers.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "PenmanMonteith" to call function
funname <- "PenmanMonteith" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "PenmanMonteith"
results <- ET(data, constants, solar="sunshine hours", wind="yes", crop="short") 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Penman-Monteith}
\keyword{evapotranspiration}
\keyword{Reference crop evapotranspiration}
\keyword{FAO-56}
\keyword{ASCE}
% __ONLY ONE__ keyword per line

