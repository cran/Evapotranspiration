\name{Evapotranspiration}
\alias{Evapotranspiration}

\title{Evapotranspiration Formulations
}
\description{
A generic function for estimating evapotranspiration from different formulations according to class of data.
}
\usage{
Evapotranspiration(data, \dots)

}
\arguments{
  \item{data}{
  A list of climate data required for estimating evapotranspiration which differs for each evapotranspiration formulations, see specific formulations for details.
  }
  \item{\dots}{
  Arguments to be passed to methods which differs for each evapotranspiration formulations, see specific formulations for details.
  }
}

\details{
  The class of data must be defined prior to call this function in order to proceed to specific formulation.
}

\author{
Danlu Guo
}

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Set data class to be "Penman" to call function
funname <- "Penman" 
class(data) <- funname 

# Call generic function Evapotranspiration(data, ...) with class "Penman"
results <- Evapotranspiration(data, constants, solar="sunshine hours", wind=
"yes", windfunction_ver = "1948", alpha = 0.08, z0 = 0.001) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{evapotranspiration}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

