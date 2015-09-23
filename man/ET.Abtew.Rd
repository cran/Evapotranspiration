\name{ET.Abtew}
\alias{ET.Abtew}

\title{Abtew Formulation
}
\description{
Implementing the Abtew formulation for estimating actual evapotranspiration.
}
\usage{
\method{ET}{Abtew}(data, constants, solar, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "Abtew" which contains the following items (climate variables) required by Abtew formulation:
  "Tmax", "Tmin", "Rs", "n" or "Cd"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Abtew formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.
  }
  \item{solar}{
  Must be either "data", "sunshine hours", "cloud" or "monthly precipitation". "data" indicates that solar radiation is to be used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{Abtew}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Abtew actual evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Abtew actual evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Abtew actual evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Abtew actual evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Abtew actual evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Abtew".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Actual Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
}
\references{
ABTEW, W. 1996. EVAPOTRANSPIRATION MEASUREMENTS AND MODELING FOR THREE WETLAND SYSTEMS IN SOUTH FLORIDA1. Wiley Online Library.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Abtew" to call function
funname <- "Abtew" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "Abtew"
results <- ET(data, constants, solar="sunshine hours")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Abtew}
\keyword{evapotranspiration}
\keyword{actual evapotranspiration}
% __ONLY ONE__ keyword per line

