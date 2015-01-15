\name{ET.JensenHaise}
\alias{ET.JensenHaise}

\title{Jensen-Haise Formulation
}
\description{
Implementing the Jensen-Haise formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{JensenHaise}(data, constants, solar, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "JensenHaise" which contains the following items (climate variables) required by Jensen-Haise formulation:
  "Tmax", "Tmin"
  }
  \item{constants}{
  A list named "constants" consists of constants required for the calculation of Penman formulation which must contain the following items:
  "Elev" - ground elevation above mean sea level in m,
  "lambda" - latent heat of vaporisation = 2.45 MJ.kg^-1,
  "lat_rad" - latitude in radians,
  "Gsc" - solar constant = 0.0820 MJ.m^-2.min^-1.
  }
  \item{solar}{
  Must be either "data", "sunshine hours", "cloud" or "monthly precipitation". "data" indicates that solar radiation data is to be used directly for calculating evapotranspiration; "sunshine hours" indicates that solar radiation is to be calculated using the real data of sunshine hours; "cloud" sunshine hours is to be estimated from cloud data; "monthly precipitation" indicates that solar radiation is to be calculated directly from monthly precipitation. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  This formulation provides a single calculation method with no alternatives available.
}

\value{
This function returns a object of class \code{JensenHaise}. Such objects are list with components:

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily Jensen-Haise potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily Jensen-Haise potential evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "Jensen-Haise".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Potential Evapotranspiration".
  }
}
\references{
JENSEN, M. E. & HAISE, H. R. 1963. Estimating evapotranspiration from solar radiation. Proceedings of the American Society of Civil Engineers, Journal of the Irrigation and Drainage Division, 89, 15-41.

}
\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "JensenHaise" to call function
funname <- "JensenHaise" 
class(data) <- funname 

# Call generic function ET(data, ...) with class "JensenHaise"
results <- ET(data, constants, solar = "sunshine hours")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Jensen-Haise}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

