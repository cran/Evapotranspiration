\name{ET.BlaneyCriddle}
\alias{ET.BlaneyCriddle}

\title{Blaney-Criddle Formulation
}
\description{
Implementing the Blaney-Criddle formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{ET}{BlaneyCriddle}(data, constants, ts="daily", solar="sunshine hours", height = F, 
message="yes", save.csv="yes", \dots) 
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Blaney-Criddle formulation:\cr
  \emph{Tmax}, \emph{Tmin} (degree Celcius), \emph{RHmin} (per cent), \emph{n} (hour) or \emph{Cd} (okta), \emph{u2} or \emph{uz} (meter per second)
  }
  \item{constants}{
    A list named \code{constants} consists of constants required for the calculation of PenPan formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{z} - height of wind instrument in m,\cr
  \emph{e0},\emph{e1},\emph{e2},\emph{e3},\emph{e4} - recommended values of 0.81917, -0.0040922, 1.0705, 0.065649, -0.0059684, -0.0005967 respectively (Table 1 in Frevert et al., 1983).
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
  }
  \item{solar}{
  Must be either \code{sunshine hours} or \code{cloud}:\cr
  \code{sunshine hours} indicates that solar radiation is to be calculated using the real data of sunshine hours;\cr
  \code{cloud} sunshine hours is to be estimated from cloud data. \cr
  Default is \code{sunshine hours}.
  }
  \item{height}{
  Must be \code{T} or \code{F}, indicating if adjustment for site elevation for arid and semi-arid regions is applied in Blaney-Criddle formulation (Allen and Brockway, 1983).\cr
  Default is \code{F} for no adjustment.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated\cr
- Evaporative surface \cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- If height adjustment has been applied on results (i.e. the value of argument \code{height})\cr
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values }
  \item{save.csv}{
  Must be either \code{yes} or \code{no}, indicating whether a .csv of ET estimates should be saved to working directory. 
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see \code{Arguments} for details.\cr
  Height adjustment for the estimations is available through argument \code{height}, please see \code{Arguments} for details.\cr
}

\value{
The function generates a list containing the following components:


  \item{ET.Daily}{
  Daily aggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  Annually ggregated estimations of Blaney-Criddle reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Blaney-Criddle reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Blaney-Criddle reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Blaney-Criddle}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Reference Crop Evapotranspiration}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message3}{
  A message to inform the users about if height adjustment has been applied to calculated Blaney-Criddle reference crop evapotranspiration.
  }
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Allen, R.G.Brockway, C.E. 1983, \emph{Estimating consumptive use on a statewide basis}. Advances in Irrigation and Drainage@ sSurviving External Pressures, ASCE, pp. 79-89.

Allen, R. & Pruitt, W. 1986. \emph{Rational Use of The FAO Blaney-Criddle Formula}. Journal of Irrigation and Drainage Engineering, 112, 139-155.

Frevert, D.K., Hill, R.W.Braaten, B.C. 1983, \emph{Estimation of FAO evapotranspiration coefficients}, Journal of Irrigation and Drainage Engineering, vol. 109, no. 2, pp. 265-270.

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

# Call ET.BlaneyCriddle under the generic function ET
results <- ET.BlaneyCriddle(data, constants, ts="daily", solar="sunshine hours", 
height= FALSE, message="yes", save.csv="yes")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{BlaneyCriddle}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

