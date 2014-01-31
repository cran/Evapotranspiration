\name{Evapotranspiration.Penpan}
\alias{Evapotranspiration.Penpan}

\title{Penpan Formulation
}
\description{
Implementing the Penpan formulation for Class-A pan evaporation.
}
\usage{
\method{Evapotranspiration}{Penpan}(data, constants, solar, alpha, overest, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "Penpan" which contains the following items (climate variables) required by Penpan formulation:
  "Tmax", "Tmin", "RHmax", "RHmin", "Rs", "n" or "Cd", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penpan formulation which must contain the following items:
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
  \item{solar}{
  Must be either "data", "sunshine hours", "cloud" or "monthly precipitation". "data" indicates that solar radiation data is to be used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface (or surrounding surface for Penpan formulation) representing the portion of the incident radiation that is reflected back at the surface. 
  }
  \item{overest}{
  Must be TRUE or FALSE, indicating if adjustment for the overestimation (i.e. divided by 1.078) of Class-A pan evaporation for Australian data is applied in Penpan formulation.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
  User-defined evaporative surface is allowed through argument \code{alpha}, please see "Arguments" for details.
  Adjustment for overestimation on the estimations are available through argument \code{height}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{Penpan}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Penpan Class-A pan evaporation.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Penpan Class-A pan evaporation.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Penpan Class-A pan evaporation.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Penpan Class-A pan evaporation.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Penpan Class-A pan evaporation.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Penpan".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Class-A Pan Evaporation".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
}
\references{

MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

ROTSTAYN, L. D., RODERICK, M. L. & FARQUHAR, G. D. 2006. A simple pan-evaporation model for analysis of climate simulations: Evaluation over Australia. Geophysical Research Letters, 33.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Penpan" to call function
funname <- "Penpan"
class(data) <- funname

# Call generic function Evapotranspiration(data, ...) with class "Penpan"
results <- Evapotranspiration(data, constants, solar="sunshine hours", alpha = 0.23, overest = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Penpan}
\keyword{evapotranspiration}
\keyword{Class-A pan evaporation}
% __ONLY ONE__ keyword per line

