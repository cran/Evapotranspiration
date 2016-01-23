\name{ET.Turc}
\alias{ET.Turc}

\title{Turc Formulation
}
\description{
Implementing the Turc formulation for estimating reference crop evapotranspiration.
}
\usage{
\method{ET}{Turc}(data, constants, ts="daily", solar="sunshine hours", 
humid=FALSE, \dots) 
}

\arguments{
  \item{data}{
  A list of data which contains the following items (climate variables) required by Turc formulation:\cr
  \emph{Tmax}, \emph{Tmin}, \emph{Rs} or \emph{n} or \emph{Cd}
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Turc formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1,\cr
  \emph{sigma} - Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.\cr
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
  \item{humid}{
  Must be \code{T} or \code{F}, indicating if adjustment for non-humid conditions is applied in Turc formulation (Alexandris et al., 2008, Equation 5b).\cr
  Default is \code{F} for no adjustment.
  }
  \item{\dots}{
  Dummy for generic function, no need to define.
  }
}

\details{
  The alternative calculation options can be selected through argument \code{solar}, please see \code{Arguments} for details.\cr
  Humidity adjustment for the estimations is available through argument \code{humid}, please see \code{Arguments} for details.\cr
}

\value{
The function prints a calculation summary to the screen containing the following elements:\cr
- ET model name and ET quantity estimated \cr
- Evaporative surface \cr
- Option for calculating solar radiation (i.e. the value of argument \code{solar})\cr
- if adjustment for non-humid conditions has been applied (i.e. the value of argument \code{humid})\cr
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values.\cr
\cr

The function also generates a list containing the following components, which is saved into a \code{csv} file named as \emph{ET_Turc.csv} in the working directory:

  \item{ET.Daily}{
  Daily aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of Turc reference crop evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily Turc reference crop evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily Turc reference crop evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{Turc}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Reference Crop Evapotranspiration}.
  }
  \item{message1}{
  A message to inform the users about how solar radiation has been calculated by using which data.
  }
  \item{message4}{
  A message to inform the users about if adjustment for non-humid conditions has been applied to calculated Turc reference crop evapotranspiration.
  }
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Turc, L. 1961, \emph{Estimation of irrigation water requirements, potential evapotranspiration: a simple climatic formula evolved up to date}. Ann. Agron, vol. 12, no. 1, pp. 13-49.

Alexandris, S., Stricevic, R.Petkovic, S. 2008, \emph{Comparative analysis of reference evapotranspiration from the surface of rainfed grass in central Serbia, calculated by six empirical methods against the Penman-Monteith formula}. European Water, vol. 21, no. 22, pp. 17-28.


}
\author{
Danlu Guo
}
\seealso{
\code{\link{ET}},\code{\link{processeddata}},\code{\link{defaultconstants}},\code{\link{constants}}
}
\examples{
# Use processed existing data set and constants from kent Town, Adelaide
data("processeddata")
data("constants")

# Call ET.Turc under the generic function ET
results <- ET.Turc(data, constants, ts="daily", 
solar="sunshine hours", humid=FALSE)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Turc}
\keyword{evapotranspiration}
\keyword{reference crop evapotranspiration}
% __ONLY ONE__ keyword per line

