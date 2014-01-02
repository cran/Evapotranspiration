\name{Evapotranspiration.Penman}
\alias{Evapotranspiration.Penman}

\title{Penman Formulation
}
\description{
Implementing the Penman formulation for estimating open-water evaporation or potential evapotranspiration
}
\usage{
\method{Evapotranspiration}{Penman}(data, constants, sunshine.hours, wind, windfunction_ver, 
alpha = 0.08, z0 = 0.001, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "Penman" which contains the following items (climate variables) required by Penman formulation:
  "Ta", Tmax", "Tmin", "n" or "Cd" or "Precip", "u2" or "uz" or "RHmean", "vabar"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penman formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1,
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "z" - height of wind instrument in m,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.
  }
  \item{sunshine.hours}{
  Must be either "data", "cloud" or "monthly precipitation". "data" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{wind}{
  Must be either "yes" or "no". "yes" indicates that the calculation will use real data of wind speed; "no" indicates that the alternative calculation without using wind data will be used in Penman formulation (Valiantzas 2006, Equation33).
  }
  \item{windfunction_ver}{
  The version of Penman wind function that will be used within the Penman formulation. Must be either "1948" or "1956". "1948" is for applying the Penman's 1948 wind function (Penman, 1948); "1956" is for applying the Penman's 1956 wind function (Penman, 1956)
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface (or surrounding surface for Penpan formulation) representing the portion of the incident radiation that is reflected back at the surface. The default is 0.08 for open-water surface which is for the calculation of Penman open-water evaporation, all other values will trigger the calculation of Penman potential evapotranspriation.
  }
  \item{z0}{
  Any value (metres), roughness height of the evaporative surface. The default is 0.08 for open-water surface which is for the calculation of Penman open-water evaporation, all other values will trigger the calculation of Penman potential evapotranspriation.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through arguments \code{sunshine.hours}, \code{wind} and \code{windfunction_ver}, please see "Arguments" for details.
  User-defined evaporative surface is allowed through arguments \code{alpha} and \code{z0}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{Penman}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Penman open-water evaporation or potential evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Penman open-water evaporation or potential evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Penman open-water evaporation or potential evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Penman open-water evaporation or potential evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Penman open-water evaporation or potential evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Penman".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is either "Open-water Evaporation" or "Potential Evaporation".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message2}{
  A message to inform the users about if actual wind data has been used in the calculations or alternative calculations has been performed without wind data, and which version of the Penman wind function has been used.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

PENMAN, H. L. 1948. Natural evaporation from open water, bare soil and grass. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences, 193, 120-145.

VALIANTZAS, J. D. 2006. Simplified versions for the Penman evaporation equation using routine weather data. Journal of Hydrology, 331, 690-702.

PENMAN, H. L. 1956. Evaporation: An introductory survey. Netherlands Journal of Agricultural Science, 4, 9-29.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Penman" to call function
funname <- "Penman" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "Penman"
results <- Evapotranspiration(data, constants, sunshine.hours="data", wind=
"yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Penman}
\keyword{evapotranspiration}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

