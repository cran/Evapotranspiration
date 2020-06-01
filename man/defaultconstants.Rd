\name{defaultconstants}
\docType{data}
\alias{defaultconstants}
\title{Universal constants Required for Calculating Evapotranspriation}
\description{
  This data set contains the universal constants required for calculating evapotranspiration in function \code{\link{ET}}, which should be kept unchanged for most conditions. Please note that additional constants may be ET models - check the manual for individual ET models for details.
}
\usage{data(defaultconstants)}
\format{A list containing 20 constant values including:
\cr
\emph{lambda} latent heat of evaporisationin = 2.45 MJ.kg^-1 at 20 degree Celcius,\cr
\emph{sigma} Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1,\cr
\emph{Gsc} solar constant = 0.0820 MJ.m^-2.min^-1\cr
\emph{Roua} mean density of air  = 1.2 kg.m^-3 at 20 degree Celcius\cr
\emph{Ca} specific heat of air = 0.001013 MJ.kg^-1.K^-1\cr
\emph{G} soil heat flux negligible for daily time-step = 0 (Allen et al., 1998, page 68)\cr
\emph{alphaA} Albedo for Class-A pan = 0.14\cr
\cr
\emph{alphaPT} Priestley-Taylor coefficient:\cr
= 1.26 for Priestley-Taylor formula (Priestley and Taylor, 1972, Sect. 6;
Eichinger et al., 1996, p.163);\cr
= 1.31 for Szilagyi-Jozsa formula (Szilagyi and Jozsa, 2008);\cr
= 1.28 for Brutsaert-Strickler formula (Brutsaert and Strickler, 1979),\cr
\cr
\emph{ap} constant in Penpan formula = 2.4,\cr
\emph{b0} constant in Morton's procedure = 1 (Chiew and McMahon, 1991, Table A1),\cr
\emph{b1} constant in Morton's procedure = 14 W.m^-2 (Chiew and McMahon, 1991, Table A1),\cr 
\emph{b2} constant in Morton's procedure = 1.2 (Chiew and McMahon, 1991, Table A1), \cr
\emph{e0} constant for Blaney-Criddle formula = 0.81917 (Frevert et al., 1983, Table 1),\cr
\emph{e1} constant for Blaney-Criddle formula = -0.0040922 (Frevert et al., 1983, Table 1),\cr
\emph{e2} constant for Blaney-Criddle formula = 1.0705 (Frevert et al., 1983, Table 1),\cr
\emph{e3} constant for Blaney-Criddle formula = 0.065649 (Frevert et al., 1983, Table 1),\cr
\emph{e4} constant for Blaney-Criddle formula = -0.0059864 (Frevert et al., 1983, Table 1),\cr
\emph{e5} constant for Blaney-Criddle formula = -0.0005967 (Frevert et al., 1983, Table 1),\cr
\emph{epsilonMo} Land surface emissivity in Morton's procedure = 0.92,\cr
\emph{sigmaMo} Stefan-Boltzmann constant in Morton's procedure = 5.67e-08 W.m^-2.K^-4.\cr
}
\seealso{
\code{\link{constants}}
}
\source{various references}
\keyword{datasets}
