\name{ET.GrangerGray}
\alias{ET.GrangerGray}

\title{Granger-Gray Formulation
}
\description{
Implementing the Granger-Gray formulation for estimating actual areal evapotranspiration.
}
\usage{
\method{ET}{GrangerGray}(data, constants, ts="daily", 
solar="sunshine hours", windfunction_ver=1948, alpha=0.23, 
message="yes", AdditionalStats="yes", save.csv="yes", \dots)
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Granger-Gray formulation:\cr
  \emph{Tmax}, \emph{Tmin} (degree Celcius), \emph{RHmax}, \emph{RHmin} (per cent), \emph{Rs} (Megajoules per sqm) or \emph{n} (hour) or \emph{Cd} (okta), \emph{u2} or \emph{uz} (meter per second)
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Granger-Gray formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1,\cr
  \emph{z} - height of wind instrument in m,\cr
  \emph{sigma} - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.\cr
  \emph{G} - soil heat flux in MJ.m^-2.day^-1, = 0 when using daily time step.
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
  \item{windfunction_ver}{
  The version of Penman wind function that will be used within the Penman formulation. Must be either \code{1948} or \code{1956}. \cr
  \code{1948} is for applying the Penman's 1948 wind function (Penman, 1948); \cr
  \code{1956} is for applying the Penman's 1956 wind function (Penman, 1956)
  Default is \code{1948}.
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface representing the portion of the incident radiation that is reflected back at the surface. \cr
  Default is 0.23 for surface covered with short reference crop.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
 - ET model name and ET quantity estimated\cr
- Evaporative surface with values of albedo\cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- The version of Penman wind function has been used (i.e. the value of argument \code{windfunction_ver})\cr
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
  The alternative calculation options can be selected through arguments \code{solar} and \code{windfunction_ver}, please see \code{Arguments} for details.\cr
  User-defined evaporative surface is allowed through argument \code{alpha}, please see \code{Arguments} for details.\cr
}

\value{
The function generates a list containing the following components:

  \item{ET.Daily}{
  Daily aggregated estimations of Granger-Gray actual areal evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Granger-Gray actual areal evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Granger-Gray actual areal evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Granger-Gray actual areal evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Granger-Gray actual areal evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Granger-Gray}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Actual Areal Evapotranspiration}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message2}{
  A message to inform the users about which version of the Penman wind function has been used.
  }
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Penman, H. L. 1948. \emph{Natural evaporation from open water, bare soil and grass}. Proceedings of the Royal Society of London. Series A. Mathematical and Physical Sciences, 193, 120-145.

Penman, H. L. 1956. \emph{Evaporation: An introductory survey}. Netherlands Journal of Agricultural Science, 4, 9-29
}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{data}},\code{\link{defaultconstants}},\code{\link{constants}},\code{\link{ET.Penman}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.GrangerGray under the generic function ET
results <- ET.GrangerGray(data, constants, ts="daily", 
solar="sunshine hours", windfunction_ver=1948, alpha=0.23, 
message="yes", AdditionalStats="yes", save.csv="yes")

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Granger-Gray}
\keyword{evapotranspiration}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

