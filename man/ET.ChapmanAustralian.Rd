\name{ET.ChapmanAustralian}
\alias{ET.ChapmanAustralian}

\title{Chapman Formulation
}
\description{
Implementing the Chapman formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{ChapmanAustralian}(data, constants, ts="daily", PenPan= T, solar="sunshine hours", alpha=0.23, \dots) 
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Chapman formulation:\cr
  \emph{Tmax}, \emph{Tmin}, \emph{RHmax}, \emph{RHmin}, \emph{Rs} or \emph{n} or \emph{Cd}, \emph{u2} or \emph{uz}

  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Chapman formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1,\cr
  \emph{z} - height of wind instrument in m,\cr
  \emph{sigma} - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,\cr
  \emph{lat} - latitude in degrees,\cr
  \emph{alphaA} - albedo for Class-A pan,\cr
  \emph{ap} - a constant in PenPan = 2.4.
  \cr
  
  The following constants are also required when argument \code{solar} has value of \code{sunshine hours}:\cr 
  \emph{as} - fraction of extraterrestrial radiation reaching earth on sunless days,\cr
  \emph{bs} - difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days.
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
  }
  \item{PenPan}{
  Must be \code{T} or \code{F}, indicating if the PenPan formulation is used for estimating Class-A pan evaporation required in Chapman formulation. If \code{T} PenPan will be used and if \code{F} the actual data of Class-A pan evaporation will be used.
  Default is \code{T} for using the PenPan formulation.
  }
  \item{solar}{
  Must be either \code{data}, \code{sunshine hours}, \code{cloud} or \code{monthly precipitation}:\cr
  \code{data} indicates that solar radiation data is to be used directly for calculating evapotranspiration; \cr
  \code{sunshine hours} indicates that solar radiation is to be calculated using the real data of sunshine hours;\cr
  \code{cloud} sunshine hours is to be estimated from cloud data; \cr
  \code{monthly precipitation} indicates that solar radiation is to be calculated directly from monthly precipitation.\cr
  Default is \code{sunshine hours}.
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of the evaporative surface incident radiation that is reflected back at the surface. \cr
  Default is 0.23 for surface covered with short reference crop.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through arguments \code{PenPan} and \code{solar}, please see \code{Arguments} for details.
}

\value{
The function prints a calculation summary to the screen containing the following elements:\cr
- ET model name and ET quantity estimated, and the value of pan coefficient (only for when potential ET is estimated)\cr
- Evaporative surface with values of albedo \cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- If the PenPan formulation is used for estimating Class-A pan evaporation required in Chapman formulation (i.e. the value of argument \code{PenPan})\cr
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values.\cr
\cr

The function also generates a list containing the following components, which is saved into a \code{csv} file named as \emph{ET_ChapmanAustralian.csv} in the working directory:

  \item{ET.Daily}{
  Daily aggregated estimations of Chapman potential evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Chapman potential evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Chapman equivalent Penmen-Monteith evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Chapman potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Chapman potential evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Chapman}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Potential Evapotranspiration}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message5}{
  A message to inform the users about if the Class-A pan evaporation is from actual data or from PenPan estimation. 
  }
}
\references{

McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Chapman, T. 2001, \emph{Estimation of evaporation in rainfall-runoff models}, in F. Ghassemi, D. Post, M. SivapalanR. Vertessy (eds), MODSIM2001: Integrating models for Natural Resources Management across Disciplines, Issues and Scales, MSSANZ, vol. 1, pp. 293-298.

}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{data}},\code{\link{defaultconstants}},\code{\link{constants}},\code{\link{ET.PenPan}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.ChapmanAustralian under the generic function ET
results <- ET.ChapmanAustralian(data, constants, ts="daily", PenPan= TRUE, 
solar="sunshine hours", alpha=0.23) 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Chapman}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

