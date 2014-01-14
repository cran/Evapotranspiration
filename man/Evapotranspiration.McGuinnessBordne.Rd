\name{Evapotranspiration.McGuinnessBordne}
\alias{Evapotranspiration.McGuinnessBordne}

\title{McGuinness-Bordne Formulation
}
\description{
Implementing the McGuinness-Bordne formulation for estimating potential evapotranspiration.
}
\usage{
\method{Evapotranspiration}{McGuinnessBordne}(data, constants, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "McGuinnessBordne" which contains the following items (climate variables) required by McGuinness-Bordne formulation:
  "Ta"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penman formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  This formulation provides a single calculation method with no alternatives available.
}

\value{
This function returns a object of class \code{McGuinnessBordne}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "McGuinness-Bordne".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Potential Evapotranspiration".
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

MCGUINNESS, J. L. & BORDNE, E. F. 1972. A comparison of lysimeter-derived potential evapotranspiration with computed values, US Dept. of Agriculture.
}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "McGuinnessBordne" to call function
funname <- "McGuinnessBordne" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "McGuinnessBordne"
results <- Evapotranspiration(data, constants)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{McGuinness-Bordne}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

