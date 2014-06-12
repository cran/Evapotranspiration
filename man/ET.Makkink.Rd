\name{ET.Makkink}
\alias{ET.Makkink}

\title{Makkink Formulation
}
\description{
Implementing the Makkink formulation for estimating potential evaporation.
}
\usage{
\method{ET}{Makkink}(data, constants, solar, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "Makkink" which contains the following items (climate variables) required by Makkink formulation:
  "Tmax", "Tmin", "Rs", "n" or "Cd"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Makkink formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days.
  }
  \item{solar}{
  Must be either "sunshine hours", "data", "cloud" or "monthly precipitation". "sunshine hours" indicates that solar radiation data will be used directly for calculating evapotranspiration; "data" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{Makkink}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Makkink potential evaporation.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Makkink potential evaporation.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Makkink potential evaporation.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Makkink potential evaporation.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Makkink potential evaporation.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Makkink".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Potential Evaporation".
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

# Set data class to be "Makkink" to call function
funname <- "Makkink" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "Makkink"
results <- ET(data, constants, solar="sunshine hours")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Makkink}
\keyword{potential evaporation}
% __ONLY ONE__ keyword per line

