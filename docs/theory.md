(theory:theory_base)=
# Underlying theory

This section contains information about the acoustic and spatial statistics theories underlying the computation implemented in Echopop. 

The [](acoustics-basics) page contains information about important acoustic quantities for deriving biological estimates in general.

The flow of computation in Echopop that combines acoustic data (in particular [NASC](NASC)) and biological data from trawls (detailed [here](apportion-trawl-stations)) to derive biomass estimates is described across a few pages:
- [](bio-estimates): derive number density from acoustic data.
- [](stratification): derive biomass density by combining. number density and animal length-weight information from trawl data.
- Kriging biomass density: perform kriging on biomass density (to be added).
- [](apportion-biomass): apportion kriged biomass density to estimate biomass for animals of different sexes, lengths, and ages.
- [](apportion-abundance): apportion abundances back-calculated from the kriged biomass estimates for animals of different sexes, lengths, and ages.

The other pages contains detailed information of spatial statistics used for kriging:
- [](semivariogram-theory): characterizing spatial variability using semivariograms.
- [](semivariogram_eq): equations of different semivariogram models.
- [](kriging-theory): interpolating observations using geostatistical estimates based on the empirical and theoretical semivariograms.