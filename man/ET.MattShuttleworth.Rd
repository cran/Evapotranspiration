\name{ET.MattShuttleworth}
\alias{ET.MattShuttleworth}

\title{Matt-Shuttleworth Formulation
}
\description{
Implementing the Matt-Shuttleworth formulation for reference crop evapotranspiration
}
\usage{
\method{ET}{MattShuttleworth}(data, constants, solar, alpha, r_s, CH, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "MattShuttleworth" which contains the following items (climate variables) required by Matt-Shuttleworth formulation:
  "Tmax", "Tmin", "RHmax", "RHmin" "n" or "Cd" or "Precip", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Matt-Schuttleworth formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "z" - height of wind instrument in m,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,
  "Roua" - mean air density = 1.20 kg.m^-3, 
  "Ca" - specific heat of air = 0.001013 MJ.kg^-1.oC^-1.
  }
  \item{solar}{
  Must be either "sunshine hours", "data", "cloud" or "monthly precipitation". "data" indicates that solar radiation data will be used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface (or surrounding surface for Penpan formulation) representing the portion of the incident radiation that is reflected back at the surface. 
  }
  \item{r_s}{
  Any value (seconds per metre), surface resistance depends on the type of reference crop.
  }
  \item{CH}{
  Any value (metres), crop height depends on the reference crop.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
  User-defined evaporative surface is allowed through arguments \code{alpha}, \code{r_s} and \code{CH}, please see "Arguments" for details.
}
\value{
This function returns a object of class \code{MattShuttleworth}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Matt-Shuttleworth".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
}
\references{

SHUTTLEWORTH, W. & WALLACE, J. 2009. Calculating the water requirements of irrigated crops in Australia using the Matt-Shuttleworth approach. Transactions of the ASABE, 52, 1895-1906.

MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "MattShuttleworth" to call function
funname <- "MattShuttleworth"
class(data) <- funname

# Call generic function ET(data, ...) with class "MattShuttleworth"
results <- ET(data, constants, solar="sunshine hours", alpha = 0.23, r_s = 70, CH = 
0.12)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Matt-Shuttleworth}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

