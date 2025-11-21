(acoustics-basics)=
# Acoustics basics

## Backscattering cross-section and target strength

### General definitions
For a given scatterer, the backscattering cross section ($\sigma_\text{bs}$, m<sup>2</sup>) is defined by {cite:p}`maclennan_2002`:

$$
    \sigma_\text{bs} = 
        \frac{r^2 I_\text{bs}(r) 10^{\alpha r/10}}{I_\text{inc}},
    \tag{1.1}
$$

where $I$ is the intensity of the incident ($\text{inc}$) and backscattered ($\text{bs}$) waves, $r$ is the range of the target (m), and $\alpha$ is the acoustic absorption coefficient (dB m<sup>-1</sup>). This quantity can be generalized using the isotropic, or spherical scattering, cross-section:

$$
    \sigma_\text{sp} = 4 \pi \sigma_\text{bs},
    \tag{1.2}
$$

where no directional dependence is assumed. For an aggregation comprising $N$ animals, the mean backscattering cross-section is:

$$
    \left< \sigma_\text{bs} \right> =
        \frac{
            \sum\limits_{j=1}^N \sigma_{\text{bs},j}
        }
        {
            N
        },
    \tag{1.3}
$$

where $\sigma_{\text{bs},j}$ is the differential backscattering cross-section of animal $j$.

This quantity can often be expressed logarithmically via target strength ($\textit{TS}$, dB re. 1 m<sup>2</sup>):

$$
    TS = 10 \log_{10} \sigma_\text{bs},
    \tag{1.4}
$$

which is the commonly used representation of individual backscattering cross-sections in the fisheries acoustics literature.

### $\textit{TS}$-length relationship

Values of $\sigma_{\text{bs},j}$ often varies as a function of the $j^\text{th}$ animal's body length, $L_j$:

$$
    \sigma_{bs,j} = \sigma_{bs,j}(L_j).
    \tag{1.5}
$$

Consequently, one common approach to estimating fish $\textit{TS}$ uses the empirical relationship between $\textit{TS}$ and $L$:

<a id="intext_eq_16"></a>

(intext_eq_16_md)=
$$
    \textit{TS} = m\log_{10}L+b,
    \tag{1.6}
$$

where $L$ is the total length, $m$ is the slope, and $b$ is the $y$-intercept. Values of $m$ and $b$ are typically species-specific; however, these relationships are often "standardized" by fixing $m = 20$ since many early experiments yielded $m$-values close to 20 {cite:p}`simmonds_maclennan_2005`. This modifies {ref}`Eq. (1.6) <eq-16>` to:

$$
    \textit{TS} = 20\log_{10}L+b_{20},
    \tag{1.7}
$$

where $b_{20}$ is the $y$-intercept associated with the fixed slope. Species-specific values for $b_{20}$ varies and are frequency-specific ({ref}`Table 1 <tbl-1>`).

(tbl-1)=
```{table} Example species-specific $\textit{TS}$-$L$ relationships.
| Species | Frequency (kHz) | $m$ | $b$ | $b_{20}$ | Ref. | Note |
|---------|-----------------|-----|---- |----------|------|------|
| *Clupea harengus* | 38 | 8.9 | -55.2 | -69.5 | {cite:t}`foote_1986` | Daytime |
| *C. harengus* | 38 | 21.2 | -74.2 | -72.5 | {cite:t}`foote_1986` | Nighttime |
| *Gadus morhua* | 38 | | | -58.8 |{cite:t}`nakken_1977` | |
| Gadoids | 38 | 18.0 | -66.2 | -68.0 | {cite:t}`foote_1987` | 
| *Merluccius productus* | 38 | | | -68.0 | {cite:t}`traynor_1996` | |
| *Pollachius virens* | 38 | | | -65.8 | {cite:t}`foote_1987` | |
| *Sprattus sprattus* | 38 | 17.2 | -60.8 | | {cite:t}`nakken_1977` | |
| *S. sprattus* | 120 | 21.4 | -55.0 | | {cite:t}`nakken_1977` | |
```

## Converting backscatter into population estimates

Methods for converting (integrated) acoustic backscatter into units of population often requires empirically measuring *in situ* $\textit{TS}$ or using models such as the regression coefficients in {ref}`Table (1) <tbl-1>`. First, the volume backscattering coefficient, $s_\text{v}$ (m<sup>-1</sup>), is computed:

$$
    s_\text{v} = 
        \frac{
            \sum\limits_{j=1}^N \sigma_{\text{bs},j}
        }
        {
            V
        },
    \tag{1.8}
$$

where $V$ is the integration volume (m<sup>3</sup>). Similar to the relationship between $\sigma_\text{bs}$ and $\textit{TS}$, $s_\text{v}$ can be expressed logarithmically via the volume backscattering strength:

$$
    S_\text{v} = 10 \log_{10} s_\text{v}.
    \tag{1.9}
$$

When $s_\text{v}$ is averaged over a finite volume of water, $S_\text{v}$ is often described as being the mean volume backscattering strength ($\textit{MVBS}$) {cite:p}`maclennan_2002`. 

Alternatively, $s_\text{v}$ can also be expressed as a function of $\left< \sigma_\text{bs} \right>$ and the volumetric animal density ($\rho_\text{v}$, animals m<sup>-3</sup>):

<a id="intext_eq_110"></a>
(intext_eq_110_md)=
$$
    s_\text{v} = \rho_\text{v} \left< \sigma_\text{bs} \right>.
    \tag{1.10}
$$

In fisheries acoustics, backscatter is often integrated over the entire water column to provide estimates per unit area instead of volume. This means that $s_\text{v}$ can be expressed as the area backscattering coefficient (m<sup>2</sup>m<sup>-2</sup>):

<a id="intext_eq_111"></a>
(intext_eq_111_md)=
$$
    s_\text{a} = 
        \int\limits_{z_1}^{z_2} s_\text{v} dz =
            s_\text{v} H,
    \tag{1.11}
$$

where $z$ is a finite depth bounded by depths $z_1$ and $z_2$, otherwise known as the integration height $H$. In practice, this vertical integration occurs over horizontal interval $x$, which is commonly referred to as an **elementary distance sampling unit (EDSU)**. This effectively breaks up continuous acoustic backscatter measurements into discrete samples (e.g. 1 nmi EDSU) that allows for easier processing and post-survey analyses (e.g. biomass estimation). This changes {ref}`Eq. (1.11) <eq-111>` to:

$$
    s_\text{a}(x) =
        \int\limits_{x_1}^{x_2} 
            \int\limits_{z_1}^{z_2}
                s_\text{v}(x, z)~ dz~ dx =
                    s_\text{v} H(x),
    \tag{1.12}
$$

where $x_1$ and $x_2$ are the along-transect distances at the start and end of the EDSU. With $s_\text{v}(x)$ integrated over the whole water column, {ref}`Eq. (1.10) <eq-110>` can be modified to compute the areal number density ($\rho_\text{a}(x)$, animals m<sup>2</sup>):

<a id="intext_eq_113"></a>
(intext_eq_113_md)=
$$
    \rho_\text{a}(x) = \frac{s_\text{a}(x)}{\left< \sigma_\text{bs} \right>}.
    \tag{1.13}
$$

It is more common to express $\rho_\text{a}$ relative to nmi<sup>2</sup> than m<sup>2</sup> or km<sup>2</sup>. This first requires converting $s_\text{a}$ into the natucal area scattering coefficient ($s_\text{A}$, or commonly abbreviated as $\textit{NASC}$) {cite:p}`foote_1994`:

$$
    s_\text{A}(x) = 4 \pi (1852)^2 s_\text{a}(x).
    \tag{1.14}
$$

Once converted, $\rho_\text{a}$ from {ref}`Eq. (1.13) <eq-113>` can be defined as:

<a id="intext_eq_115"></a>
(intext_eq_115_md)=
$$
    \rho_\text{A} =
        \frac{s_\text{A}}{4 \pi \left< \sigma_\text{bs} \right>} =
            \frac{s_\text{A}}{\sigma_\text{sp}}
    \tag{1.15}
$$

where $\rho_\text{A}(x)$ is in units of animals nmi<sup>-2</sup>.