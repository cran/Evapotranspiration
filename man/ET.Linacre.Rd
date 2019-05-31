\name{ET.Linacre}
\alias{ET.Linacre}

\title{Linacre Formulation
}
\description{
Implementing the Linacre formulation for estimating actual evapotranspiration.
}
\usage{
\method{ET}{Linacre}(data, constants, ts="daily", message="yes", AdditionalStats="yes", 
save.csv="yes", \dots) 
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Linacre formulation:\cr
  \emph{Tmax}, \emph{Tmin}, \emph{Tdew} (degree Celcius)
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Linacre formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lat} - latitude in degrees.\cr
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
  \item{AdditionalStats}{
  "yes" or "no" indicating whether monthly averaged and annual averaged ET should be calculated.
  }
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
  Daily aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Linacre actual evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Linacre actual evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Linacre actual evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Linacre}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Actual Evapotranspiration}.
  }
}
\references{
Linacre, E. T. 1977. \emph{A simple formula for estimating evaporation rates in various climates, using temperature data alone}. Agricultural meteorology, 18, 409-424.

}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{data}},\code{\link{defaultconstants}},\code{\link{constants}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.Linacre under the generic function ET
results <- ET.Linacre(data, constants, ts="daily",
message="yes", AdditionalStats="yes", save.csv="yes") 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Linacre}
\keyword{evapotranspiration}
\keyword{actual evapotranspiration}
% __ONLY ONE__ keyword per line

