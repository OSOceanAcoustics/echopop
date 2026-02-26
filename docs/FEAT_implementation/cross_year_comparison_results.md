(comparison_results)=

# Cross-year comparisons

The figure below summarizes the signed percent differences in total abundance, biomass, and $S_\text{A}$ between EchoPro and Echopop across all survey years, for both the transect and kriged population estimates. Each cell represents:

$$
    \text{signed percent difference} = \text{sign}(\Delta) \cdot \frac{|\Delta|}{\frac{\text{EchoPro} + \text{Echopop}}{2}} \times 100,
$$

where $\Delta = \text{EchoPro} - \text{Echopop}$. Positive values (red) indicate EchoPro exceeds Echopop; negative values (blue) indicate the opposite. Cells close to zero (white) indicate strong agreement between the two implementations.

![Cross-year comparisons](../_static/comparisons/cross_year_comparisons_20260225.png)

:::{admonition} What do these outputs represent?
:class: note
Importantly, the values shown in the figure above correspond to the outputs of Echopop and EchoPro where the following conditions were enforced:
- **Age-2+ only**
- **Biodata processed using KS strata**
- **Kriged estimates were extrapolated over the entire mesh grid**
:::

## Results

### Transect estimates

Agreement between EchoPro and Echopop transect estimates is generally strong across all survey years, with differences becoming increasingly small in more recent surveys.

#### 1995 - 2009, 2012 - 2021
Percent differences in abundance, biomass, and $S_\text{A}$ are all below 1.0%, indicating near-identical results between the two implementations for this period. The small observed differences can likely be attributed to the [nuanced changes to how Echopop](#echopro-vs-echopop) computes the number and weight proportions (e.g., unsexed fish are not assumed to be female) compared to EchoPro.

#### 2011
Differences remain small but are slightly more variable. In 2011, $S_\text{A}$ and biomass show strong agreement while abundance is marginally elevated at 1.5%. This can be attributed to three hauls (51, 52, 53) being treated as Canadian in the EchoPro biodata, and as American in the dataset Echopop ingests. In the latter, the Canadian haul offset (100) is not added to these particular hauls, so they are not included in the biodata. This is because hauls 151, 152, and 153 are contained within the associated stratification files, **but not** 51, 52, and 53.

When these hauls are re-incorporated as Canadian hauls, the difference in abundance (1.5%) and biomass (0.3%) both decrease to values less than 0.1%.

### Kriged estimates

#### 1995, 2015-2021
Percent differences in kriged abundance, biomass, and $S_\text{A}$ are all below 1.0%, indicating near-identical results between the two implementations for this period. The small observed differences can likely be attributed to the [nuanced changes to how Echopop](#echopro-vs-echopop) handles the kriging algorithm compared to EchoPro.

#### 1998, 2001, 2005 - 2013

#### 2003