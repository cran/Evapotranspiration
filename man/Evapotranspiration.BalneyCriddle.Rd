\name{Evapotranspiration.BalneyCriddle}
\alias{Evapotranspiration.BalneyCriddle}

\title{Balney-Criddle Formulation
}
\description{
Implementing the Balney-Criddle formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{Evapotranspiration}{BalneyCriddle}(data, constants, sunshine.hours, height, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "BalneyCriddle" which contains the following items (climate variables) required by Balney-Criddle formulation:
  "Ta", "RHmin", "n" or "Cd", "u2" or "uz"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Balney-Criddle formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lat_rad" - latitude in radians,
  "z" - height of wind instrument in m,
  "e0", "e1", "e2", "e3", "e4" - recommended values of 0.81917, -0.0040922, 1.0705, 0.065649, -0.0059684, -0.0005967 respectively (Table 1 in Frevert et al., 1983).
  }
  \item{sunshine.hours}{
  Must be either "data" or "cloud". "data" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data.
  }
  \item{height}{
  Must be TRUE or FALSE, indicating if adjustment for site elevation for arid and semi-arid regions is applied in Balney-Criddle formulation (Allen and Brockway, 1983).
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{sunshine.hours}, please see "Arguments" for details.
  Height adjustment for the estimations is available through argument \code{height}, please see "Arguments" for details.
}

\value{
This function returns a object of class \code{BalneyCriddle}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Balney-Criddle reference crop evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Balney-Criddle reference crop evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Balney-Criddle reference crop evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Balney-Criddle reference crop evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Balney-Criddle reference crop evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Balney-Criddle".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message3}{
  A message to inform the users about if height adjustment has been applied to calculated Balney-Criddle reference crop evapotranspiration.
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

ALLEN, R. G. & BROCKWAY, C. E. Estimating consumptive use on a statewide basis.  Advances in Irrigation and Drainage@ sSurviving External Pressures, 1983. ASCE, 79-89.

FREVERT, D. K., HILL, R. W. & BRAATEN, B. C. 1983. Estimation of FAO evapotranspiration coefficients. Journal of irrigation and drainage engineering, 109, 265-270.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "BalneyCriddle" to call function
funname <- "BalneyCriddle" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "BalneyCriddle"
results <- Evapotranspiration(data, constants, sunshine.hours="data", height = TRUE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BalneyCriddle}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

