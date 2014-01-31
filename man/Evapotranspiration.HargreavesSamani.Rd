\name{Evapotranspiration.HargreavesSamani}
\alias{Evapotranspiration.HargreavesSamani}

\title{Hargreaves-Samani Formulation
}
\description{
Implementing the Hargreaves-Samani formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{Evapotranspiration}{HargreavesSamani}(data, constants, \dots)
}

\arguments{
  \item{data}{
  A list of data in class "HargreavesSamani" which contains the following items (climate variables) required by Hargreaves-Samani formulation:
  "Tmax", "Tmin"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Hargreaves-Samani formulation which must contain the following items:
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
This function returns a object of class \code{HargreavesSamani}. Such objects are list with components:

  \item{PET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{PET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{PET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{PET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{PET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{PET_formulation}{
  A character string containing the name of the formulation used which equals to "Hargreaves-Samani".
  }
  \item{PET_type}{
  A character string containing the type of the estimation obtained which is "Reference Crop Evapotranspiration".
  }
}
\references{
MCMAHON, T., PEEL, M., LOWE, L., SRIKANTHAN, R. & MCVICAR, T. 2012. Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

HARGREAVES, G. H. & SAMANI, Z. A. 1985. Reference crop evapotranspiration from ambient air temperature. American Society of Agricultural Engineers.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "HargreavesSamani" to call function
funname <- "HargreavesSamani" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "HargreavesSamani"
results <- Evapotranspiration(data, constants)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hargreaves-Samani}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

