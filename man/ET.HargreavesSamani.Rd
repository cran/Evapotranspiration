\name{ET.HargreavesSamani}
\alias{ET.HargreavesSamani}

\title{Hargreaves-Samani Formulation
}
\description{
Implementing the Hargreaves-Samani formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{ET}{HargreavesSamani}(data, constants, ts="daily",
message="yes", AdditionalStats="yes", save.csv="yes", \dots)
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Hargreaves-Samani formulation:
  \emph{Tmax}, \emph{Tmin} (degree Celcius)
  }
  \item{constants}{
    A list named \code{constants} consists of constants required for the calculation of Hargreaves-Samani formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1.\cr
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated \cr
- Evaporative surface with values of albedo \cr
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
  Daily aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Hargreaves-Samani reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Hargreaves-Samani}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Reference Crop Evapotranspiration}.
  }
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Hargreaves, G.H.Samani, Z.A. 1985, \emph{Reference crop evapotranspiration from ambient air temperature}. American Society of Agricultural Engineers.
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

# Call ET.HargreavesSamani under the generic function ET
results <- ET.HargreavesSamani(data, constants, ts="daily", message="yes", 
AdditionalStats="yes", save.csv="yes")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Hargreaves-Samani}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

