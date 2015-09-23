\name{ET.Linacre}
\alias{ET.Linacre}

\title{Linacre Formulation
}
\description{
Implementing the Linacre formulation for estimating actual evapotranspiration.
}
\usage{
\method{ET}{Linacre}(data, constants, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "Linacre" which contains the following items (climate variables) required by Linacre formulation:
  "Tmax", "Tmin" and "Tdew"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Linacre formulation which must contain the following items:
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
This function returns a object of class \code{Linacre}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Linacre actual evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Linacre actual evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Linacre".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Actual Evapotranspiration".
  }
}
\references{
LINACRE, E. T. 1977. A simple formula for estimating evaporation rates in various climates, using temperature data alone. Agricultural meteorology, 18, 409-424.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Linacre" to call function
funname <- "Linacre" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "Linacre"
results <- ET(data, constants)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linacre}
\keyword{evapotranspiration}
\keyword{actual evapotranspiration}
% __ONLY ONE__ keyword per line

