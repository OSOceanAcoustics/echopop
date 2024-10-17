(echopro-vs-echopop)=
# Differences between the `EchoPro` and `Echopop`

There are documented differences between the `EchoPro` MATLAB program and `Echopop` workflows. Changes include new features (or changes to methods in `EchoPro`) not previously implemented in `Echopop` (<span style="font-size:5mm;">✨</span>), bugfixes (<span style="font-size:5mm;">🐛</span>), language-specific differences between MATLAB and Python (<span style="font-size:5mm;">👽</span>), features `EchoPro` that have not yet been implemented in `Echopop` (<span style="font-size:5mm;">🧩</span>), and features in `EchoPro` that are not currently in development for `Echopop` (<span style="font-size:5mm;">📕</span>).

## Data ingestion
- <span style="font-size:6mm;">[✨]</span> The creation of external `*.xlsx` transect-region-haul key files from Echoview `*.xlsx` exports are no longer required
- <span style="font-size:6mm;">[🧩]</span> Manually mapping transect numbers, regions, and haul numbers together from Echoview `*.xlsx` exports is currently not implemented

## Transect analysis
- <span style="font-size:6mm;">[🐛]</span> Age bins can now be manually defined in the configuration YAML file and do not automatically subsume out-of-bounds values (e.g. $\alpha >= 21$ years all being included within $\alpha = 20$)
- <span style="font-size:6mm;">[🐛]</span> Sex-specific TS-length regressions no longer arbitrarily redefine 'unsexed' fish as being 'female'
- <span style="font-size:6mm;">[👽]</span> Interpolated weights across $\ell$ can produce different results in `EchoPro` and `Echopop` due to differences in numerical precision and how missing data are handled
- <span style="font-size:6mm;">[🧩]</span> Incorporating net selectivity into biological distributions (across $\ell$ and $\alpha$) is currently not implemented
- <span style="font-size:6mm;">[✨/🧩]</span> Filtering out specific haul numbers, ages, sexes, and other plausible contrasts have not yet been implemented
- <span style="font-size:6mm;">[📕]</span> Support for incorporating observer data within `Echopop` is currently not supported

## Stratified analysis
- <span style="font-size:6mm;">[✨]</span> Confidence intervals at user-defined significance levels are now provided for CV. This includes bootstrap/resampling bias calculations and several methods for computing the confidence intervals: empirical, percentile, standard, $t$-standard, $t$-jackknife, bias-corrected (BC), and bias-corrected and accelerated (BCa)
- <span style="font-size:6mm;">[✨]</span> Estimators and corresponding confidence intervals are now provided for (resampled) mean abundance/biomass
- <span style="font-size:6mm;">[🧩]</span> The Jolly and Hampton (1990) resampling approach for computing overall variability in $\textit{NASC}$ estimates has not been implemented in `Echopop`. It also currently limited to being computed only over INPFC strata

## General spatial methods
- <span style="font-size:6mm;">[✨]</span> A new alternative method for cropping the kriging mesh has been implemented based on the convex hull shape of the survey transects  (i.e. `crop_method="convex_hull"`)
- <span style="font-size:6mm;">[✨/🐛]</span> The implementation of the `EchoPro` kriged mesh cropping method relied on interpolating the transect line extents of three discrete regions to account for the island of Haida Gwaii. In `EchoPro`, these are manually defined. `Echopop` defaults to this implementation (i.e. `crop_method="transect_ends"`) by discretizing the transect lines into each of these three regions based on their respective headings in cardinal directions 
- <span style="font-size:6mm;">[🐛]</span> Modifications were made to prevent erroneous gaps in the survey coverage shape that subsequently included unexpected mesh nodes


## Variogram analysis
- <span style="font-size:6mm;">[✨]</span> The semivariogram fitting GUI (i.e. `Survey.variogram_gui()`) has been enabled to allow for user changes to optimization algorithm parameters (e.g. `gradient_tolerance`)
- <span style="font-size:6mm;">[👽]</span> Differences between the non-linear least squares optimization algorithm in MATLAB and Python (via `lmfit`) differs in a number of ways that can produce different parameter estimates. These differences can range from precision error to much larger; however, the relative importance placed on parameter order remains the same in `Echopop` as it does in `EchoPro`
  
## Kriging analysis
- <span style="font-size:6mm;">[✨]</span> Kriged unaged biomass estimates are reapportioned along length and age based on the distributions computed from aged fish. When certain length-bins are missing, `EchoPro` imputes values using the closest length bin corresponding to either male or female fish. `Echopop` does this imputation using the closest sex-specific length bins instead
- <span style="font-size:6mm;">[🧩]</span> `Echopop` currently does not support kriged abundance or $\textit{NASC}$ back-calculation from biomass estimates
- <span style="font-size:6mm;">[🧩]</span> Kriging has only been fully tested, validated, and implemented for biomass