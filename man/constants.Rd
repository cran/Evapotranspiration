\name{constants}
\docType{data}
\alias{constants}
\title{Constants Required for Calculating Evapotranspriation}
\description{
  This data set contains the constants required for calculating evapotranspiration in function \code{ET} for Kent Town station in Adelaide, Australia.
}
\usage{climatedata}
\format{A list containing 46 constant values including:
"lambda", latent heat of evaporisationin = 2.45 MJ.kg^-1 at 20 degree Celcius,
"sigma", Stefan-Boltzmann constant = 4.903*10^-9 MJ.K^-4.m^-2.day^-1.
"Gsc", solar constant = 0.0820 MJ.m^-2.min^-1,
"lat", latitude = -34.9211 degrees for Kent Town station,
"lat_rad", latitude in radians = -0.6095 radians for Kent Town station,
"as", fraction of extraterrestrial radiation reaching earth on sunless days = 0.23 for Australia (Roderick, 1999, page 181),
"bs", difference between fracion of extraterrestrial radiation reaching full-sun days and that on sunless days = =0.5 for Australia (Roderick, 1999, page 181),
"Elev", ground elevation above mean sea level = 48m for Kent Town station,
"z", height of wind instrument = 10m for Kent Town station,
"Roua", mean density of air	= 1.2 kg.m^-3 at 20 degree Celcius,
"Ca", specific heat of air = 0.001013 MJ.kg^-1.K^-1,
"G", soil heat flux negligible for daily time-step = 0 (Allen et al., 1998, page 68),
"alphaA", Albedo for Class-A pan = 0.14,
"alphaPT", Priestley-Taylor coefficient = 1.26 for Priestley-Taylor formula (Priestley and Taylor, 1972, Sect. 6;
Eichinger et al., 1996, p.163); = 1.31 for Szilagyi-Jozsa formula (Szilagyi and Jozsa, 2008); = 1.28 for Brutsaert-Strickler formula (Brutsaert and Strickler, 1979),
"ap", constant in Penpan formula = 2.4,
"fz", constant in Morton's procedure = 28.0 W.m^-2.mbar^-1 for CRAE model for T >= 0 degree Celcius, = 28.0*1.15 W.m^-2.mbar^-1 for CRAE model for T < 0 degree Celcius; = 25.0 W.m^-2.mbar^-1 for CRWE model for T >= 0 degree Celcius, = 28.75 W.m^-2.mbar^-1 for CRWE model for T < 0 degree Celcius (Morton, 1983a, page65),
"b0", constant in Morton's procedure = 1 (Chiew and McMahon, 1991, Table A1),
"b1", constant in Morton's procedure = 14 W.m^-2 (Chiew and McMahon, 1991, Table A1), 
"b2", constant in Morton's procedure = 1.2 (Chiew and McMahon, 1991, Table A1), 
"a_0", constant for estimating sunshine hours from cloud cover data = 11.9 for Adelaide (Chiew and McMahon, 1991, Table A1),
"b_0", constant for estimating sunshine hours from cloud cover data = -0.15 for Adelaide, 
"c_0", constant for estimating sunshine hours from cloud cover data = -0.25 for Adelaide,
"d_0", constant for estimating sunshine hours from cloud cover data = -0.0107 for Adelaide
"e0", constant for Blaney-Criddle formula = 0.81917 (Frevert et al., 1983, Table 1),
"e1", constant for Blaney-Criddle formula = -0.0040922 (Frevert et al., 1983, Table 1),
"e2", constant for Blaney-Criddle formula = 1.0705 (Frevert et al., 1983, Table 1),
"e3", constant for Blaney-Criddle formula = 0.065649 (Frevert et al., 1983, Table 1),
"e4", constant for Blaney-Criddle formula = -0.0059864 (Frevert et al., 1983, Table 1),
"e5", constant for Blaney-Criddle formula = -0.0005967 (Frevert et al., 1983, Table 1),
"gammaps", product of Psychrometric constant and atmospheric pressure as sea level = 0.66 mbar. degree Celcius^-1 for CRAE model for T >= 0 degree Celcius, = 0.66/1.15 mbar. degree Celcius^-1 for CRAE model for T < 0 degree Celcius,
"epsilonMo", Land surface emissivity in Morton's procedure = 0.92,
"PA", annual precipitation = 285.8mm for Kent Town station, 
"alphaMo", constant in Morton's procedure = 17.27 when T >= 0 degree Celcius, = 21.88 when T < 0 degree Celcius,
"betaMo", constant in Morton's procedure = 237.3 degree Celcius when T >= 0 degree Celcius, = 265.5 degree Celcius when T < 0 degree Celcius,
"sigmaMo", Stefan-Boltzmann constant in Morton's procedure = 5.67e-08 W.m^-2.K^-4,
"lambdaMo", latent heat of vaporisation in Morton's procedure = 28.5W.day.kg^-1 when T >= 0 degree Celcius, = 28.5*1.15W.day.kg^-1 when T < 0 degree Celcius.
}
\source{various references}
\keyword{datasets}
