\name{ET.MattShuttleworth}
\alias{ET.MattShuttleworth}

\title{Matt-Shuttleworth Formulation
}
\description{
Implementing the Matt-Shuttleworth formulation for reference crop evapotranspiration
}
\usage{
\method{ET}{MattShuttleworth}(data, constants, ts="daily", solar="sunshine hours", 
alpha=0.23, r_s=70, CH=0.12, message="yes", AdditionalStats="yes", save.csv="yes", \dots)
}

\arguments{
  \item{data}{
  A list which contains the following items (climate variables) required by Matt-Shuttleworth formulation:\cr
  \emph{Tmax}, \emph{Tmin} (degree Celcius), \emph{RHmax}, \emph{RHmin} (per cent), \emph{Rs} (Megajoules per sqm) or \emph{n} (hour) or \emph{Cd} (okta), \emph{u2} or \emph{uz} (meter per second)
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Matt-Shuttleworth formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1,\cr
  \emph{z} - height of wind instrument in m,\cr
  \emph{sigma} - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,\cr
  \emph{Roua} - mean air density = 1.20 kg.m^-3,\cr
  \emph{Ca} - specific heat of air = 0.001013 MJ.kg^-1.oC^-1.\cr
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
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface representing the portion of the incident radiation that is reflected back at the surface. \cr
  Default is 0.23 for surface covered with short reference crop, which is for the calculation of Matt-Shuttleworth reference crop evaporation.
  }
  \item{r_s}{
  Any value (seconds per metre), surface resistance depends on the type of reference crop.\cr
  Default is 70 for short reference crop.
  }
  \item{CH}{
  Any value (metres), crop height depends on the reference crop.\cr
  Default is 0.12 for short reference crop.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated\cr
- Evaporative surface with values of albedo, surface resistance and crop height\cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
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
  The alternative calculation options can be selected through argument \code{solar}, please see \code{Arguments} for details.\cr
  User-defined evaporative surface is allowed through arguments \code{alpha}, \code{r_s} and \code{CH}, please see \code{Arguments} for details.
}
\value{
The function generates a list containing the following components:
  \item{ET.Daily}{
  Daily aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Matt-Shuttleworth reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \emph{Matt-Shuttleworth}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \emph{Reference Crop Evapotranspiration}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
}
\references{

Shuttleworth, W. & Wallace, J. 2009. \emph{Calculating the water requirements of irrigated crops in Australia using the Matt-Shuttleworth approach}. Transactions of the ASABE, 52, 1895-1906.

McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

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

# Call ET.MattShuttleworth under the generic function ET
results <- ET.MattShuttleworth(data, constants, ts="daily",
solar="sunshine hours", alpha=0.23, r_s=70, CH=0.12,
message="yes", AdditionalStats="yes", save.csv="yes")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Matt-Shuttleworth}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

