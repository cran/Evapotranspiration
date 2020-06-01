\name{ET}
\alias{ET}

\title{ET Formulations
}
\description{
A generic function including 17 different specific methods that are all named following the format of \emph{ET.methodname}. Once specific function is called the corresponding calculations are performed and a calculation summary is printed to screen.

}
\usage{
ET(data, constants, \dots)

}
\arguments{
  \item{data}{
  A list of climate data required for estimating evapotranspiration which differs for each evapotranspiration formulations, see specific formulations for details.
  }
  \item{constants}{
  A list named \code{constants} consists of constants required for the ET models which mdiffer for specific ET models - refer to the manual for individual models for details.}
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

\examples{
# Use processed existing data set from kent Town, Adelaide
data("processeddata")
data("constants")

# Call generic function ET() - leads to the use of Penman model
results_default <- ET(processeddata, constants, save.csv="no")

# Call generic function ET() - leads to the use of Penman-Monteith model
results_crop <- ET(processeddata, constants, crop = "short", save.csv="no")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{evapotranspiration}
\keyword{open-water evaporation}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line
