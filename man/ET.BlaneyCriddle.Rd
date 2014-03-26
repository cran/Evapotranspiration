\name{ET.BlaneyCriddle}
\alias{ET.BlaneyCriddle}

\title{Blaney-Criddle Formulation
}
\description{
Implementing the Blaney-Criddle formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{ET}{BlaneyCriddle}(data, constants, solar, height, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "BlaneyCriddle" which contains the following items (climate variables) required by Blaney-Criddle formulation:
  "Tmax", "Tmin", "RHmin", "n", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Blaney-Criddle formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lat_rad" - latitude in radians,
  "z" - height of wind instrument in m,
  "e0", "e1", "e2", "e3", "e4" - recommended values of 0.81917, -0.0040922, 1.0705, 0.065649, -0.0059684, -0.0005967 respectively (Table 1 in Frevert et al., 1983).
  }
  \item{solar}{
  Must be either "sunshine hours" or "cloud". "sunshine hours" indicates that solar radiation is to be calculated using the data of sunshine hours; for "cloud" sunshine hours is to be estimated from cloud data.
  }
  \item{height}{
  Must be TRUE or FALSE, indicating if adjustment for site elevation for arid and semi-arid regions is applied in Blaney-Criddle formulation (Allen and Brockway, 1983).
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see "Arguments" for details.
  Height adjustment for the estimations is available through argument \code{height}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{BlaneyCriddle}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Blaney-Criddle reference crop evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Blaney-Criddle reference crop evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Blaney-Criddle".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message3}{
  A message to inform the users about if height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

ALLEN, R. & PRUITT, W. 1986. Rational Use of The FAO Blaney-Criddle Formula. Journal of Irrigation and Drainage Engineering, 112, 139-155.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "BlaneyCriddle" to call function
funname <- "BlaneyCriddle" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "BlaneyCriddle"
results <- ET(data, constants, solar="sunshine.hours", height = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BlaneyCriddle}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

