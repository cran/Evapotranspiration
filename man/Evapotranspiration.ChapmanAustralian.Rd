\name{Evapotranspiration.ChapmanAustralian}
\alias{Evapotranspiration.ChapmanAustralian}

\title{Chapman Formulation
}
\description{
Implementing the Chapman formulation for estimating equivalent Penmen-Monteith reference crop evapotranspiration.
}
\usage{
\method{Evapotranspiration}{ChapmanAustralian}(data, constants, Penpan, solar, alpha, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "ChapmanAustralian" which contains the following items (climate variables) required by Chapman formulation:
  "Epan" or "Tmax" and "Tmin" and "RHmax" and "RHmin" and "u2" or "uz" and "Rs" or "n" or "Cd"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penman formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "z" - height of wind instrument in m,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,
  "lat" - latitude in degrees,
  "alphaA" - albedo for Class-A pan,
  "ap" - a constant in Penpan = 2.4.
  }
  \item{Penpan}{
  Must be TRUE or FALSE, indicating if the Penpan formulation is used for estimating Class-A pan evaporation required in Chapman formulation. If TRUE Penpan will be used and if FALSE the actual data of Class-A pan evaporation will be used.
  }
  \item{solar}{
  Must be either "data", "sunshine hours", or "cloud". "data" indicates that solar radiation data is used directly for calculating evapotranspiration, "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data.
  }
  \item{alpha}{
  Must be defined if Penpan = TRUE. Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface (or surrounding surface for Penpan formulation) representing the portion of the incident radiation that is reflected back at the surface. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through arguments \code{Penpan} and \code{solar}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{ChapmanAustralian}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Chapman equivalent Penmen-Monteith reference crop evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Chapman equivalent Penmen-Monteith reference crop evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Chapman equivalent Penmen-Monteith reference crop evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Chapman equivalent Penmen-Monteith reference crop evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Chapman equivalent Penmen-Monteith reference crop evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Chapman".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Equivalent Penmen-Monteith Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message5}{
  A message to inform the users about if the Class-A pan evaporation is from actual data or from Penpan estimation. 
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "ChapmanAustralian" to call function
funname <- "ChapmanAustralian" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "ChapmanAustralian"
results <- Evapotranspiration(data, constants, Penpan = TRUE, solar = "sunshine hours", 
alpha = 0.23)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Chapman}
\keyword{evapotranspiration}
\keyword{equivalent Penmen-Monteith reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

