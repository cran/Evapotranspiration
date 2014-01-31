\name{Evapotranspiration.BrutsaertStrickler}
\alias{Evapotranspiration.BrutsaertStrickler}

\title{Brutsaert-Strickler Formulation
}
\description{
Implementing the Brutsaert-Strickler formulation for actual areal evapotranspiration
}
\usage{
\method{Evapotranspiration}{BrutsaertStrickler}(data, constants, solar, alpha, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "BrutsaertStrickler" which contains the following items (climate variables) required by Brutsaert-Strickler formulation:
  "Tmax", "Tmin", "RHmax", "RHmin", "RS", "n" or "Cd", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Brutsaert-Strickler formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "z" - height of wind instrument in m,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.
  }
  \item{solar}{
  Must be either "data", "sunshine hours", cloud" or "monthly precipitation". "data" indicates that solar radiation data is used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface (or surrounding surface for Penpan formulation) representing the portion of the incident radiation that is reflected back at the surface. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
  User-defined evaporative surface is allowed through argument \code{alpha}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{BrutsaertStrickler}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Brutsaert-Strickler actual areal evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Brutsaert-Strickler actual areal evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Brutsaert-Strickler actual areal evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Brutsaert-Strickler actual areal evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Brutsaert-Strickler actual areal evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Brutsaert-Strickler".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Actual Areal Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
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

# Set data class to be "BrutsaertStrickler" to call function
funname <- "BrutsaertStrickler"
class(data) <- funname

# Call generic function Evapotranspiration(data, ...) with class "BrutsaertStrickler"
results <- Evapotranspiration(data, constants, solar="sunshine.hours", alpha = 0.08)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Brutsaert-Strickler}
\keyword{evapotranspiration}
\keyword{actual areal evapotranspiration}
% __ONLY ONE__ keyword per line

