\name{ET.Turc}
\alias{ET.Turc}

\title{Turc Formulation
}
\description{
Implementing the Turc formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{ET}{Turc}(data, constants, solar, humid, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "Turc" which contains the following items (climate variables) required by Turc formulation:
  "Tmax", "Tmin", "Rs", "n" or "Cd"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Turc formulation which must contain the following items:
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
  \item{humid}{
  Must be TRUE or FALSE, indicating if adjustment for non-humid conditions is applied in Turc formulation (Alexandris et al., 2008, Equation 5b).
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
  Humidity adjustment for the estimations is available through argument \code{humid}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{Turc}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Turc reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Turc reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Turc".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message4}{
  A message to inform the users about if adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

TURC, L. 1961. Estimation of irrigation water requirements, potential evapotranspiration: a simple climatic formula evolved up to date. Ann. Agron, 12, 13-49.

ALEXANDRIS, S., STRICEVIC, R. & PETKOVIC, S. 2008. Comparative analysis of reference evapotranspiration from the surface of rainfed grass in central Serbia, calculated by six empirical methods against the Penman-Monteith formula. Eur Water, 21, 17-28.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Turc" to call function
funname <- "Turc" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "Turc"
results <- ET(data, constants, solar="sunshine hours", humid = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Turc}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

