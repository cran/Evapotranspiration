\name{ET.JensenHaise}
\alias{ET.JensenHaise}

\title{Jensen-Haise Formulation
}
\description{
Implementing the Jensen-Haise formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{JensenHaise}(data, constants, ts="daily", solar="sunshine hours", \dots)
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Jensen-Haise formulation:
  \emph{Tmax}, \emph{Tmin}, \emph{Rs} or \emph{n} or \emph{Cd}
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Jensen-Haise formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1.\cr
  \cr
  
  The following constants are also required when argument \code{solar} has value of \code{sunshine hours}:\cr 
  \emph{as} - fraction of extraterrestrial radiation reaching earth on sunless days,\cr
  \emph{bs} - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days.
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
  }
  \item{solar}{
  Must be either \code{data}, \code{sunshine hours}, \code{cloud} or \code{monthly precipitation}:\cr
  \code{data} indicates that solar radiation data is to be used directly for calculating evapotranspiration; \cr
  \code{sunshine hours} indicates that solar radiation is to be calculated using the real data of sunshine hours;\cr
  \code{cloud} sunshine hours is to be estimated from cloud data; \cr
  \code{monthly precipitation} indicates that solar radiation is to be calculated directly from monthly precipitation.\cr
  Default is \code{sunshine hours}.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  This formulation provides a single calculation method with no alternatives available.
}

\value{
The function prints a calculation summary to the screen containing the following elements:\cr
- ET model name and ET quantity estimated \cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values.\cr
\cr

The function also generates a list containing the following components, which is saved into a \code{csv} file named as \emph{ET_JensenHaise.csv} in the working directory:

  \item{ET.Daily}{
  Daily aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Jensen-Haise potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Jensen-Haise potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Jensen-Haise potential evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Jensen-Haise}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Potential Evapotranspiration}.
  }
}
\references{
Jensen, M.E.Haise, H.R. 1963, \emph{Estimating evapotranspiration from solar radiation}. Proceedings of the American Society of Civil Engineers, Journal of the Irrigation and Drainage Division, vol. 89, pp. 15-41.

Prudhomme, C.Williamson, J. 2013, \emph{Derivation of RCM-driven potential evapotranspiration for hydrological climate change impact analysis in Great Britain: a comparison of methods and associated uncertainty in future projections}. Hydrol. Earth Syst. Sci., vol. 17, no. 4, pp. 1365-1377.

Xu, C.Y.Singh, V.P. 2000, \emph{Evaluation and generalization of radiation-based methods for calculating evaporation}., Hydrological Processes, vol. 14, no. 2, pp. 339-349.

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

# Call ET.JensenHaise under the generic function ET
results <- ET.JensenHaise(data, constants, ts="daily", solar="sunshine hours")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Jensen-Haise}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

