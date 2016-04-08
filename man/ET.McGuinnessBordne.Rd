\name{ET.McGuinnessBordne}
\alias{ET.McGuinnessBordne}

\title{McGuinness-Bordne Formulation
}
\description{
Implementing the McGuinness-Bordne formulation for estimating potential evapotranspiration.
}
\usage{
\method{ET}{McGuinnessBordne}(data, constants, ts="daily", \dots) 
}

\arguments{
  \item{data}{
    A list of data which contains the following items (climate variables) required by McGuinness-Bordne formulation:
  \emph{Tmax}, \emph{Tmin}
  }
  \item{constants}{
    A list named \code{constants} consists of constants required for the calculation of Jensen-Haise formulation which must contain the following items:\cr
  \emph{Elev} - ground elevation above mean sea level in m,\cr
  \emph{lambda} - latent heat of vaporisation = 2.45 MJ.kg^-1,\cr
  \emph{lat_rad} - latitude in radians,\cr
  \emph{Gsc} - solar constant = 0.0820 MJ.m^-2.min^-1.\cr
  }
  \item{ts}{
  Must be either \code{daily}, \code{monthly} or \code{annual}, which indicates the disired time step that the output ET estimates should be on. 
  Default is \code{daily}.
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
- Time step of the output ET estimates (i.e. the value of argument \code{ts})\cr
- Units of the output ET estimates\cr
- Time duration of the ET estimation\cr
- Number of ET estimates obtained in the entire time-series\cr
- Basic statistics of the estimated ET time-series including \emph{mean}, \emph{max} and \emph{min} values.\cr
\cr

The function also generates a list containing the following components, which is saved into a \code{csv} file named as \emph{ET_McGuinnessBordne.csv} in the working directory:

  \item{ET.Daily}{
  Daily aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.Monthly}{
  Monthly aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.Annual}{
  Annually aggregated estimations of McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.MonthlyAve}{
  Monthly averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET.AnnualAve}{
  Annually averaged estimations of daily McGuinness-Bordne potential evapotranspiration.
  }
  \item{ET_formulation}{
  Name of the formulation used which equals to \code{McGuinness-Bordne}.
  }
  \item{ET_type}{
  Type of the estimation obtained which is \code{Potential Evapotranspiration}.
  }
}
\references{
Oudin, L., Hervieu, F., Michel, C., Perrin, C., Andreassian, V., Anctil, F.Loumagne, C. 2005, \emph{Which potential evapotranspiration input for a lumped rainfall-runoff model?: Part 2-Towards a simple and efficient potential evapotranspiration model for rainfall-runoff modelling}. Journal of Hydrology, vol. 303, no. 1-4, pp. 290-306.

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

# Call ET.McGuinnessBordne under the generic function ET
results <- ET.McGuinnessBordne(data, constants, ts="daily")
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{McGuinness-Bordne}
\keyword{evapotranspiration}
\keyword{potential evapotranspiration}
% __ONLY ONE__ keyword per line

