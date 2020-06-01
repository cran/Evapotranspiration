\name{ET.MortonCRWE}
\alias{ET.MortonCRWE}

\title{Morton CRWE Formulation
}
\description{
Implementing the Morton CRWE formulation for estimating potential evapotranspiration or shallow lake evaporation.
}
\usage{
\method{ET}{MortonCRWE}(data, constants, ts="monthly", est="potential ET",
solar="sunshine hours", Tdew= T, alpha = NULL, message="yes", AdditionalStats="yes",
save.csv="no", \dots)
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Morton CRWE formulation:\cr
  \emph{Tmax}, \emph{Tmin}, \emph{Tdew} (degree Celcius) or \emph{va} or \emph{RHmax} and \emph{RHmin}, \emph{Rs} (Megajoules per sqm) or \emph{n} (hour) or \emph{Cd} (okta)
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Morton CRWE formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{PA} - annual precipitation in mm, required when precipitation data is not available,\cr
  \emph{sigma} - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,\cr
  \emph{lat} - latitude in degrees,\cr
 \emph{epsilonMo} - surface emissivity = 0.92 (Morton, 1986),\cr
  \emph{fz} - A constant in Morton's procedure = 25.0 Wm^-2.mbar^-1 for T >= 0 degree Celcius, and = 28.75 Wm^-2.mbar^-1 for T >= 0 degree Celcius for CRWE model (Morton, 1986),\cr
  \emph{b0} - A constants in Morton's procedure, = 1.12 for CRWE model, (Morton, 1986)
  \emph{b1} - A constant in Morton's procedure, = 13 for CRWE model (Morton, 1986),\cr
  \emph{b2} - A constant in Morton's procedure, = 1.12 for CRWE model (Morton, 1986),\cr
  \emph{gammaps} - Produce of Psychrometric constant and atmospheric pressure as sea level, = 0.66 mbar. degree Celcius^-1 for T >= 0 degree Celcius, = 0.66/1.15 mbar. degree Celcius^-1 for T < 0 degree Celcius (Morton, 1983),\cr
  \cr
  \emph{alphaMo} - a constant in Morton's procedure, = 17.27 when T >= 0 degree Celcius, = 21.88 when T < 0 degree Celcius (Morton, 1983),\cr
  \emph{betaMo} - a constant in Morton's procedure, = 237.3 degree Celcius  when T >= 0 degree Celcius, = 265.5 degree Celcius, when T < 0 degree Celcius (Morton, 1983),\cr
  \emph{sigmaMo} - Stefan-Boltzmann constant in Morton's procedure, = 5.67e-08 W.m^-2.K^-4 (Morton, 1983),\cr
  \emph{lambdaMo} - Latent heat of vaporisation in Morton's procedure, = 28.5W.day.kg^-1 when T >= 0 degree Celcius, = 28.5*1.15W.day.kg^-1 when T < 0 degree Celcius,\cr
  }
  \item{ts}{
  Must be either \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on.
  Default is \code{monthly}.
  }
  \item{solar}{
  Must be either \code{data}, \code{sunshine hours}, \code{cloud} or \code{monthly precipitation}:\cr
  \code{data} indicates that solar radiation data is to be used directly for calculating evapotranspiration; \cr
  \code{sunshine hours} indicates that solar radiation is to be calculated using the real data of sunshine hours;\cr
  \code{cloud} sunshine hours is to be estimated from cloud data; \cr
  \code{monthly precipitation} indicates that solar radiation is to be calculated directly from monthly precipitation.\cr
  Default is \code{sunshine hours}.
  }
  \item{est}{
  Must be either \code{potential ET} or \code{shallow lake ET}:\cr
  \code{potential ET} proceeds to estimating potential evapotranspiration;\cr
  \code{shallow lake ET} proceeds to estimating shallow lake evaporation.\cr
  Default is \code{potential ET}.
  }
  \item{Tdew}{
   Must be \code{T} or \code{F}, indicating if real data of dew point temperature is used for calculating the radiation in Morton's formulations, if \code{T} the data will be used and if \code{F} the dew point temperature will be calculated from data of daily vapour pressure.
   Default is \code{T} for using actual dew point temperature data.
  }
  \item{alpha}{
  Only needed if argument \code{solar} has value of \code{data}. \cr
 Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface representing the portion of the incident radiation that is reflected back at the surface. \cr
  Default is \code{NULL} in line with the default use of sunshine hours to estimate solar radiation (i.e. argument \code{solar} is \code{sunshine hours}.
  }
  \item{message}{
  Must be either \code{yes} or \code{no}, indicating whether message should be printed for calculation summary including the following elements:\cr
- ET model name and ET quantity estimated (i.e. the value of argument \code{est}) \cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- If the actual dew point temperature data are used (i.e. the value of argument \code{Tdew})\cr
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
  The type of evapotranspiration calculated can be selected through argument \code{est}, please see \code{Arguments} for details.
  The alternative calculation options can be selected through argument \code{solar} and \code{Tdew}, please see \code{Arguments} for details.
}

\value{
The function generates a list containing the following components:

  \item{ET.Daily}{
  Daily aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily MortonCRWE potential evapotranspiration or shallow lake evaporation.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{MortonCRWE}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is either \code{Potential Evapotranspiration} or \code{Shallow Lake Evaporation}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message6}{
  A message to inform the users about if actual dew point temperature has been used in the calculations or alternative calculations has been performed without dew point temperature data.
  }
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Morton, F.I. 1983, \emph{Operational estimates of lake evaporation}. Journal of Hydrology, vol. 66, no. 1-4, pp. 77-100.
}
\author{
Danlu Guo
}
\seealso{
\code{\link{processeddata}},\code{\link{defaultconstants}},\code{\link{constants}},\code{\link{ET.MortonCRWE}}
}
\examples{
# Use processed existing data set and constants from
# kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.MortonCRWE under the generic function ET
results <- ET.MortonCRWE(processeddata, constants, ts="monthly",
est="potential ET", solar="sunshine hours", Tdew= TRUE,
alpha = NULL, message="yes", AdditionalStats="yes", save.csv="no")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{MortonCRWE}
\keyword{shallow lake}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

