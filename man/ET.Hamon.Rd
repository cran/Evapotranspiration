\name{ET.Hamon}
\alias{ET.Hamon}

\title{Hamon Formulation
}
\description{
Implementing the Hamon formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{Hamon}(data, constants = NULL, ts="daily", message="yes", save.csv="yes", \dots) 
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Hamon formulation:\cr
  \emph{Tmax}, \emph{Tmin} (degree Celcius), \emph{n} (hour)
  }
  \item{constants}{
  Dummy argument with a \code{NULL} value.
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated\cr
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values.}
  \item{save.csv}{
  Must be either \code{yes} or \code{no}, indicating whether a .csv of ET estimates should be saved to working directory. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  This formulation provides a single calculation method with no alternatives available.
}

\value{
The function generates a list containing the following components:

  \item{ET.Daily}{
  Daily aggregated estimations of Hamon potential evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Hamon potential evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Hamon potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Hamon potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Hamon potential evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Hamon}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Potential Evapotranspiration}.
  }
}
\references{
Hamon, W. R. 1961. \emph{Estimating potential evapotranspiration}. Journal of the Hydraulics Division, 87, 107-120.

Oudin, L., Hervieu, F., Michel, C., Perrin, C., Andreassian, V., Anctil, F.Loumagne, C. 2005, \emph{Which potential evapotranspiration input for a lumped rainfall-runoff model?: Part 2-Towards a simple and efficient potential evapotranspiration model for rainfall-runoff modelling}. Journal of Hydrology, vol. 303, no. 1-4, pp. 290-306.


}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{data}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.Hamon under the generic function ET
results <- ET.Hamon(data, ts="daily", message="yes", save.csv="yes") 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hamon}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

