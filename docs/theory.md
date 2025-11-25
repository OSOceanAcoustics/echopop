(theory:theory_base)=
# Underlying theory

This section contains information about the acoustic and spatial statistics theories underlying the computation implemented in Echopop. 

The [](acoustics-basics) page contains information about important acoustic quantities for deriving biological estimates in general.

The flow of computation in Echopop that combines acoustic data and biological data from trawls to derive biomass estimates is described across a few pages:
- [](bio-estimates): derive number density from acoustic data.
- [](stratification): derive biomass density by combining. number density and animal length-weight information from trawl data.

The other pages contains detailed information of spatial statistics used for kriging:
- [](semivariogram-theory): characterizing spatial variability using semivariograms.
- [](semivariogram_eq): equations of different semivariogram models.
- [](kriging-theory): interpolating observations using geostatistical estimates based on the empirical and theoretical semivariograms.