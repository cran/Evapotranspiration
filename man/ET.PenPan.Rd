\name{ET.PenPan}
\alias{ET.PenPan}

\title{PenPan Formulation
}
\description{
Implementing the PenPan formulation for Class-A pan evaporation.
}
\usage{
\method{ET}{PenPan}(data, constants, ts="daily", solar="sunshine hours",
alpha=0.23, est="potential ET", pan_coeff=0.71, overest= F, message="yes",
AdditionalStats="yes", save.csv="no", \dots)
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by PenPan formulation:\cr
  \emph{Tmax}, \emph{Tmin} (degree Celcius), \emph{RHmax}, \emph{RHmin} (per cent), \emph{Rs} (Megajoules per sqm) or \emph{n} (hour) or \emph{Cd} (okta), \emph{u2} or \emph{uz} (meter per second)
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of PenPan formulation which must contain the following items:\cr
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
  \item{solar}{
  Must be either \code{data}, \code{sunshine hours}, \code{cloud} or \code{monthly precipitation}:\cr
  \code{data} indicates that solar radiation data is to be used directly for calculating evapotranspiration; \cr
  \code{sunshine hours} indicates that solar radiation is to be calculated using the real data of sunshine hours;\cr
  \code{cloud} sunshine hours is to be estimated from cloud data; \cr
  \code{monthly precipitation} indicates that solar radiation is to be calculated directly from monthly precipitation.\cr
  Default is \code{sunshine hours}.
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of surface surrounding the evaporation pan representing the portion of the incident radiation that is reflected back at the surface. \cr
  Default is 0.23 for surface covered with short reference crop.
  }
  \item{overest}{
  Must be \code{T} or \code{F}, indicating if adjustment for the overestimation (i.e. divided by 1.078) of Class-A pan evaporation for Australian data is applied in PenPan formulation.\cr
  Default is \code{F} for no adjustment.
  }
  \item{est}{
  Must be either \code{pan} or \code{potential ET} to specify if estimation for the Class-A pan evaporation or potential evapotranspriation is performed.\cr
  Default is \code{potential ET} for estimating potential evapotranspriation.
  }
  \item{pan_coeff}{
  Only required if argument \code{est} has value of \code{potential ET}, which defines the pan coefficient used to adjust the estimated pan evaporation to the potential ET required.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated (i.e. the value of argument \code{est}), and the value of pan coefficient (only for when potential ET is estimated)\cr
- Evaporative surface with values of albedo \cr
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
  User-defined evaporative surface is allowed through argument \code{alpha}, please see \code{Arguments} for details.\cr
  Adjustment for overestimation on the estimations are available through argument \code{height}, please see \code{Arguments} for details.\cr
}

\value{
The function generates a list containing the following components:

  \item{ET.Daily}{
  Daily aggregated estimations of PenPan Class-A pan evaporation/potential evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of PenPan Class-A pan evaporation/potential evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of PenPan Class-A pan evaporation/potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily PenPan Class-A pan evaporation/potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily PenPan Class-A pan evaporation/potential evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{PenPan}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Class-A Pan Evaporation} or \code{Potential Evapotranspiration} depending on the value of \code{est}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
}
\references{

McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Rotstayn, L. D., Roderick, M. L. & Farquhar, G. D. 2006. \emph{A simple pan-evaporation model for analysis of climate simulations: Evaluation over Australia}. Geophysical Research Letters, 33.
}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{processeddata}},\code{\link{defaultconstants}},\code{\link{constants}}
}
\examples{
# Use processed existing data set and constants from
# kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.PenPan under the generic function ET
results <- ET.PenPan(processeddata, constants, ts="daily",
solar="sunshine hours", alpha=0.23,
est="potential ET", pan_coeff=0.71, overest= FALSE,
message="yes", AdditionalStats="yes", save.csv="no")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{PenPan}
\keyword{evapotranspiration}
\keyword{Class-A pan evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line
