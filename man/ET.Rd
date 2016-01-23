\name{ET}
\alias{ET}

\title{ET Formulations
}
\description{
A generic function including 17 different specific methods that are all named following the format of \emph{ET.methodname}. Once specific function is called the corresponding calculations are performed and a calculation summary is printed to screen.

}
\usage{
ET(data, constants, ts="daily", crop=NULL, alpha=NULL, \dots)

}
\arguments{
  \item{data}{
  A list of climate data required for estimating evapotranspiration which differs for each evapotranspiration formulations, see specific formulations for details.
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the calculation of Penman-Monteith formulation which must contain the following items:\cr
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
  \item{crop}{
  Only for the Penman-Monteith model, must be either \code{short} or \code{tall}.\cr 
  \code{short} indicates that the method for FAO-56 hypothetical short grass will be applied (Allen et al., 1998, Equation 6); \cr
  \code{tall} indicates that the method for  ASCE-EWRI Standardised crop will be applied (ASCE, 2005, Equation 1, Table 1).\cr
  Default is \code{NULL}.
  }
  \item{alpha}{
  Any numeric value between 0 and 1 (dimensionless), albedo of evaporative surface representing the portion of the incident radiation that is reflected back at the surface. \cr
  Default is \code{NULL}.
  }
  \item{\dots}{
  Arguments to be passed to methods which differs for each evapotranspiration formulations, see specific formulations for details.
  }
}

\details{
  Individual ET methods can be called by substituting the 'methodname' by the function name (e.g. \code{\link{ET.Penman}} to call the Penman model).\cr
  \cr
  When the ET model selection is not specified by users, this function determines the default model to use based on the availability of climate data presented. Wherever data are available, the more comprehensive, physically-based models are always preferred over the empirical models, in the following hierarchy:\cr
  \cr
- If all variables of \emph{Tmax/Tmin} and \emph{RHmax/RHmin} and either \emph{uz} or \emph{u2}, and either \emph{Rs} of \emph{n} or \emph{Cd} are available, and short crop surface is specified in argument: \cr 
  Penman-Monteith FAO56 (\code{\link{ET.PenmanMonteith}} with \emph{crop = "short"});\cr
  \cr 
- If all variables of \emph{Tmax/Tmin} and \emph{RHmax/RHmin} and either \emph{uz} or \emph{u2}, and either \emph{Rs} of \emph{n} or \emph{Cd} are available, and long crop surface is specified in argument: \cr 
  Penman-Monteith ASCE-EWRI (\code{\link{ET.PenmanMonteith}} with \emph{crop = "long"});\cr
  \cr 
- If all variables of \emph{Tmax/Tmin} and \emph{RHmax/RHmin} and either \emph{uz} or \emph{u2}, and either \emph{Rs}      of \emph{n} or \emph{Cd} are available, and no surface is specified:\cr 
  Penman (\code{\link{ET.Penman}});\cr
  \cr 
- If all variables of \emph{Tmax/Tmin} and \emph{RHmax/RHmin}, and either \emph{Rs} of \emph{n} or \emph{Cd} are available:\cr 
  Priestley-Taylor (\code{\link{ET.PriestleyTaylor}});\cr
  \cr 
- If all variables of \emph{Tmax/Tmin} and either \emph{Rs} of \emph{n} or \emph{Cd} are available:\cr 
  Makkink (\code{\link{ET.Makkink}});\cr
  \cr 
- If all variables of \emph{Tmax/Tmin} are available:\cr 
  Hargreaves-Samani (\code{\link{ET.HargreavesSamani}}).
  
}

\author{
Danlu Guo
}
\seealso{
\code{\link{ET.PenmanMonteith}},\code{\link{ET.Penman}},\code{\link{ET.PriestleyTaylor}},\code{\link{ET.Makkink}},\code{\link{ET.HargreavesSamani}},\code{\link{processeddata}},\code{\link{defaultconstants}},\code{\link{constants}}
}
\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Call generic function ET without specifying 
# specific method
results <- ET(data, constants, ts="daily")
# triggers the use of Penman model according to the data availability

results <- ET(data, constants, ts="daily", crop="short")
# triggers the use of Penman-Monteith model according to the data availability
# evaporative surface

# if no wind data present
data$uz <- NULL
results <- ET(data, constants, ts="daily", alpha=0.23)
# triggers the use of Priestley-Taylor model according to the data availability
# evaporative surface

# if no wind data & RH data present
data$uz <- data$RHmax <- data$RHmin <- NULL
results <- ET(data, constants, ts="daily")
# triggers the use of Makkink model according to the data availability

# if only temperature data present
data$uz <- data$RHmax <- data$RHmin <- data$n <- NULL
results <- ET(data, constants, ts="daily")
# triggers the use of Hargreaves-Samani model according to the data availability
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{evapotranspiration}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line
