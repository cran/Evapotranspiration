\name{ET.Romanenko}
\alias{ET.Romanenko}

\title{Romanenko Formulation
}
\description{
Implementing the Romanenko formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{Romanenko}(data, constants, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "Romanenko" which contains the following items (climate variables) required by Romanenko formulation:
  "Tmax", "Tmin", "RHmax" and "RHmin"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Romanenko formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1
  "as" - fraction of extraterrestrial radiation reaching earth on sunless days,
  "bs" - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days,
  "sigma" - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  This formulation provides a single calculation method with no alternatives available.
}

\value{
This function returns a object of class \code{Romanenko}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Romanenko potential evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Romanenko potential evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Romanenko potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Romanenko potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Romanenko potential evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Romanenko".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Potential Evapotranspiration".
  }
}
\references{
OUDIN, L., HERVIEU, F., MICHEL, C., PERRIN, C., ANDREASSIAN, V., ANCTIL, F. & LOUMAGNE, C. 2005. Which potential evapotranspiration input for a lumped rainfall-runoff model?: Part 2-Towards a simple and efficient potential evapotranspiration model for rainfall-runoff modelling. Journal of Hydrology, 303, 290-306.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Romanenko" to call function
funname <- "Romanenko" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "Romanenko"
results <- ET(data, constants)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Romanenko}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

