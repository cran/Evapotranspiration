\name{constants}
\alias{constants}
\docType{data}
\title{Constants Required for Calculating Evapotranspriation}
\description{
   This data set contains the universal constants, and examples of other variable constants required for calculating evapotranspiration in function \code{\link{ET}}, based on the climatic condition at Kent Town station in Adelaide, Australia.
}
\usage{data(constants)}
\format{A list containing 36 constant values including:\cr
\cr
- 20 universal constants, which should be kept unchanged for most conditions:\cr
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
*Note: a re-calibrated value of 13.4 W.m^-2 was recommended to achieve achieve a Priestley-Taylor coefficient of 1.26 (Wang et al., 2009), rather the original value (14 W.m^-2) used by Morton that gave a Priestley-Taylor coefficient of 1.32;\cr

\emph{b2} constant in Morton's procedure = 1.2 (Chiew and McMahon, 1991, Table A1), \cr
*Note: a re-calibrated value of 1.13 was recommended to achieve achieve a Priestley-Taylor coefficient of 1.26 (Wang et al., 2009), rather the original value (1.2) used by Morton that gave a Priestley-Taylor coefficient of 1.32;\cr

\emph{e0} constant for Blaney-Criddle formula = 0.81917 (Frevert et al., 1983, Table 1),\cr
\emph{e1} constant for Blaney-Criddle formula = -0.0040922 (Frevert et al., 1983, Table 1),\cr
\emph{e2} constant for Blaney-Criddle formula = 1.0705 (Frevert et al., 1983, Table 1),\cr
\emph{e3} constant for Blaney-Criddle formula = 0.065649 (Frevert et al., 1983, Table 1),\cr
\emph{e4} constant for Blaney-Criddle formula = -0.0059864 (Frevert et al., 1983, Table 1),\cr
\emph{e5} constant for Blaney-Criddle formula = -0.0005967 (Frevert et al., 1983, Table 1),\cr
\emph{epsilonMo} Land surface emissivity in Morton's procedure = 0.92,\cr
\emph{sigmaMo} Stefan-Boltzmann constant in Morton's procedure = 5.67e-08 W.m^-2.K^-4.\cr
\cr
- 16 variable constants, which are specific for the climatic condition at Kent Town station in Adelaide, Australia:\cr
\cr
\emph{lat} latitude = -34.9211 degrees for Kent Town station,\cr
\emph{lat_rad} latitude in radians = -0.6095 radians for Kent Town station,\cr
\emph{as} fraction of extraterrestrial radiation reaching earth on sunless days = 0.23 for Australia (Roderick, 1999, page 181),\cr
\emph{bs} difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days = 0.5 for Australia (Roderick, 1999, page 181),\cr
\emph{Elev} ground elevation above mean sea level = 48m for Kent Town station,\cr
\emph{z} height of wind instrument = 10m for Kent Town station,\cr
\cr
\emph{fz} constant in Morton's procedure: \cr
= 28.0 W.m^-2.mbar^-1 for CRAE model for T >= 0 degree Celcius;\cr
*Note: a re-calibrated value of 29.2 W.m^-2.mbar^-1 was recommended to achieve achieve a Priestley-Taylor coefficient of 1.26 (Wang et al., 2009), rather the original value (28.0 W.m^-2.mbar^-1) used by Morton that gave a Priestley-Taylor coefficient of 1.32;\cr

= 28.0*1.15 W.m^-2.mbar^-1 for CRAE model for T < 0 degree Celcius; \cr

= 25.0 W.m^-2.mbar^-1 for CRWE model for T >= 0 degree Celcius;\cr
= 28.75 W.m^-2.mbar^-1 for CRWE model for T < 0 degree Celcius (Morton, 1983a, page65).\cr
\cr
\emph{a_0} constant for estimating sunshine hours from cloud cover data = 11.9 for Adelaide (Chiew and McMahon, 1991, Table A1),\cr
\emph{b_0} constant for estimating sunshine hours from cloud cover data = -0.15 for Adelaide, \cr
\emph{c_0} constant for estimating sunshine hours from cloud cover data = -0.25 for Adelaide,\cr
\emph{d_0} constant for estimating sunshine hours from cloud cover data = -0.0107 for Adelaide,
\emph{gammaps} product of Psychrometric constant and atmospheric pressure as sea level:\cr
= 0.66 mbar. degree Celcius^-1 for CRAE model for T >= 0 degree Celcius;\cr
= 0.66/1.15 mbar. degree Celcius^-1 for CRAE model for T < 0 degree Celcius. \cr
\emph{PA} annual precipitation = 285.8mm for Kent Town station, \cr
\cr
\emph{alphaMo} constant in Morton's procedure:\cr 
= 17.27 when T >= 0 degree Celcius;\cr
= 21.88 when T < 0 degree Celcius.\cr

\emph{betaMo} constant in Morton's procedure:\cr
= 237.3 degree Celcius when T >= 0 degree Celcius;\cr
= 265.5 degree Celcius when T < 0 degree Celcius.\cr

\emph{lambdaMo} latent heat of vaporisation in Morton's procedure:\cr
= 28.5W.day.kg^-1 when T >= 0 degree Celcius;\cr
= 28.5*1.15W.day.kg^-1 when T < 0 degree Celcius.
}
\seealso{
\code{\link{defaultconstants}}
}
\references{
McMahon, T., Peel, M., Lowe, L., Srikanthan, R. & McVicar, T. 2012. \emph{Estimating actual, potential, reference crop and pan evaporation using standard meteorological data: a pragmatic synthesis}. Hydrology and Earth System Sciences Discussions, 9, 11829-11910.

Allen, R. G., Pereira, L. S., Raes, D. & Smith, M. 1998. \emph{Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage}. paper 56. FAO, Rome, 300, 6541.

Szilagyi, J., & Jozsa, J. 2008. New findings about the complementary relationship-based evaporation estimation methods. Journal of Hydrology, 354(1-4), 171-186.

Brutsaert, W., & Stricker, H. 1979. An advection-aridity approach to estimate actual regional evapotranspiration. Water Resources Research, 15(2), 443-450. 

Chiew, F. H. S., & McMahon, T. A. 1991. The applicability of Morton's and Penman's evapotranspiration estimates in rainfall-runoff modelling. JAWRA Journal of the American Water Resources Association, 27(4), 611-620. 

Frevert, D.K., Hill, R.W.Braaten, B.C. 1983, Estimation of FAO evapotranspiration coefficients, Journal of Irrigation and Drainage Engineering, vol. 109, no. 2, pp. 265-270.

Roderick, M. L. 1999. Estimating the diffuse component from daily and monthly measurements of global radiation. Agricultural and Forest Meteorology, 95(3), 169-185. 

Wang, Q. J., McConachy, F. L. N., Chiew, F. H. S., James, R., de Hoedt, G. C., & Wright, W. J. 2009. Maps of Evapotranspiration. Retrieved from Melbourne, Australia: http://www.bom.gov.au/climate/averages/climatology/evapotrans/text/et-description.pdf

Morton, F. I. 1983. Operational estimates of areal evapotranspiration and their significance to the science and practice of hydrology. Journal of Hydrology, 66(1-4), 1-76. doi:http://dx.doi.org/10.1016/0022-1694(83)90177-4
}

\keyword{datasets}
