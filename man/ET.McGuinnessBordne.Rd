\name{ET.McGuinnessBordne}
\alias{ET.McGuinnessBordne}

\title{McGuinness-Bordne Formulation
}
\description{
Implementing the McGuinness-Bordne formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{McGuinnessBordne}(data, constants, \dots) 
}

\arguments{
  \item{data}{
  A list of data in class "McGuinnessBordne" which contains the following items (climate variables) required by McGuinness-Bordne formulation:
  "Tmax", "Tmin"
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

  \item{ET.Daily}{
  A \code{zoo} object containing daily aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.Monthly}{
  A \code{zoo} object containing monthly aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.Annual}{
  A \code{zoo} object containing annually aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  A \code{zoo} object containing monthly averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  A \code{zoo} object containing annually averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET_formulation}{
  A character string containing the name of the formulation used which equals to "McGuinness-Bordne".
  }
  \item{ET_type}{
  A character string containing the type of the estimation obtained which is "Potential Evapotranspiration".
  }
}
\references{
XU, C. Y. & SINGH, V. P. 2000. Evaluation and generalization of radiation-based methods for calculating evaporation. Hydrological Processes, 14, 339-349.
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

# Call generic function ET(data, ...) with class "McGuinnessBordne"
results <- ET(data, constants)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{McGuinness-Bordne}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

