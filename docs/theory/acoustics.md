(acoustics-basics)=
# Acoustics basics

## Backscattering cross-section and target strength

### General definitions
For a given scatterer, the backscattering cross section ($\sigma_\text{bs}$, m<sup>2</sup>) is defined by[^maclennan_et_al]:

$$
\sigma_\text{bs} = 
\frac{r^2 I_\text{bs}(r) 10^{\alpha r/10}}{I_\text{inc}},
\tag{1}
$$

where $I$ is the intensity of the incident ($\text{inc}$) and backscattered ($\text{bs}$) waves, $r$ is the range of the target (m), and $\alpha$ is the acoustic absorption coefficient (dB m<sup>-1</sup>). This quantity can be generalized using the isotropic, or spherical scattering, cross-section:

$$
\sigma_\text{sp} = 4 \pi \sigma_\text{bs},
\tag{2}
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
\tag{3}
$$

where $\sigma_{\text{bs},j}$ is the differential backscattering cross-section of animal $j$.

This quantity can often be expressed logarithmically via target strength ($\textit{TS}$, dB re. 1 m<sup>2</sup>):

$$
TS = 10 \log_{10} \sigma_\text{bs},
\tag{4}
$$

which is the commonly used representation of individual backscattering cross-sections in the fisheries acoustics literature.

### $\textit{TS}$-length relationship

Values of $\sigma_{\text{bs},j}$ often varies as a function of the $j^\text{th}$ animal's body length, $L_j$:

$$
\sigma_{bs,j} = \sigma_{bs,j}(L_j).
\tag{5}
$$

Consequently, one common approach to estimating fish $\textit{TS}$ uses the empirical relationship between $\textit{TS}$ and $L$:

$$
\textit{TS} = m\log_{10}L+b,
\tag{6},
\label{eq:TS-L}
$$

where $L$ is the total length, $m$ is the slope, and $b$ is the $y$-intercept. Values of $m$ and $b$ are typically species-specific; however, these relationships are often "standardized" by fixing $m = 20$ since many early experiments yielded $m$-values close to 20[^fisheries_acoustics]. This modifies $\eqref{eq:TS-L}$ to:

$$
\textit{TS} = 20\log_{10}L+b_{20},
\tag{7}
$$

where $b_{20}$ is the $y$-intercept associated with the fixed slope. Species-specific values for $b_{20}$ varies and are frequency-specific ({numref}`tbl:ts_length_params`).

```{table} Example species-specific $\textit{TS}$-$L$ relationships.
:name: tbl:ts_length_params
| Species | Frequency (kHz) | $m$ | $b$ | $b_{20}$ | Ref. | Note |
|---------|-----------------|-----|---- |----------|------|------|
| *Clupea harengus* | 38 | 8.9 | -55.2 | -69.5 | [^foote_1986] | Daytime |
| *C. harengus* | 38 | 21.2 | -74.2 | -72.5 | [^foote_1986] | Nighttime |
| *Gadus morhua* | 38 | | | -58.8 | [^nakken_olsen] | |
| Gadoids | 38 | 18.0 | -66.2 | -68.0 | [^foote_1987] | 
| *Merluccius productus* | 38 | | | -68.0 | [^traynor] | |
| *Pollachius virens* | 38 | | | -65.8 | [^foote_1987] | |
| *Sprattus sprattus* | 38 | 17.2 | -60.8 | | [^nakken_olsen] | |
| *S. sprattus* | 120 | 21.4 | -55.0 | | [^nakken_olsen] | |
```

## Converting backscatter into population estimates

Methods for converting (integrated) acoustic backscatter into units of population often requires empirically measuring *in situ* $\textit{TS}$ or using models such as the regression coefficients in {numref}`tbl:ts_length_params`. First, the volume backscattering coefficient, $s_\text{v}$ (m<sup>-1</sup>), is computed:

$$
s_\text{v} = 
\frac{
    \sum\limits_{j=1}^N \sigma_{\text{bs},j}
}
{
    V
},
\tag{8}
$$

where $V$ is the integration volume (m<sup>3</sup>). Similar to the relationship between $\sigma_\text{bs}$ and $\textit{TS}$, $s_\text{v}$ can be expressed logarithmically via the volume backscattering strength:

$$
S_\text{v} = 10 \log_{10} s_\text{v}.
\tag{9}
$$

When $s_\text{v}$ is averaged over a finite volume of water, $S_\text{v}$ is often described as being the mean volume backscattering strength ($\textit{MVBS}$)[^maclennan_et_al]. 

Alternatively, $s_\text{v}$ can also be expressed as a function of $\left< \sigma_\text{bs} \right>$ and the volumetric animal density ($\rho_\text{v}$, animals m<sup>-3</sup>):

$$
s_\text{v} = \rho_\text{v} \left< \sigma_\text{bs} \right>.
\tag{10}
\label{eq:sv_pv_relationship}
$$

In fisheries acoustics, backscatter is often integrated over the entire water column to provide estimates per unit area instead of volume. This means that $s_\text{v}$ can be expressed as the area backscattering coefficient (m<sup>2</sup>m<sup>-2</sup>):

$$
s_\text{a} = 
\int\limits_{z_1}^{z_2} s_\text{v} dz =
s_\text{v} H,
\tag{11}
\label{eq:sa_integrate}
$$

$$ s_\text{a} \approx \bar{s}_\text{v} H $$

where $z$ is a finite depth bounded by depths $z_1$ and $z_2$, otherwise known as the integration height $H$. In practice, this vertical integration occurs over horizontal interval $x$, which is commonly referred to as an **elementary distance sampling unit (EDSU)**. This effectively breaks up continuous acoustic backscatter measurements into discrete samples (e.g. 1 nmi EDSU) that allows for easier processing and post-survey analyses (e.g. biomass estimation). This changes $\eqref{eq:sa_integrate}$ to:

$$
s_\text{a}(x) =
\int\limits_{x_1}^{x_2} 
\int\limits_{z_1}^{z_2}
s_\text{v}(x, z)~ dz~ dx =
s_\text{v} H(x),
\tag{12}
$$

where $x_1$ and $x_2$ are the along-transect distances at the start and end of the EDSU. With $s_\text{v}(x)$ integrated over the whole water column, $\eqref{eq:sv_pv_relationship}$ can be modified to compute the areal number density ($\rho_\text{a}(x)$, animals m<sup>2</sup>):

$$
\rho_\text{a}(x) = \frac{s_\text{a}(x)}{\left< \sigma_\text{bs} \right>}.
\tag{13}
\label{eq:areal_density_m}
$$

It is more common to express $\rho_\text{a}$ relative to nmi<sup>2</sup> than m<sup>2</sup> or km<sup>2</sup>. This first requires converting $s_\text{a}$ into the natucal area scattering coefficient ($s_\text{A}$, or commonly abbreviated as $\textit{NASC}$)[^foote_knudsen_1994]:

$$
s_\text{A}(x) = 4 \pi (1852)^2 s_\text{a}(x).
\tag{14}
$$

Once converted, $\rho_\text{a}$ from $\eqref{eq:areal_density_m}$ can be defined as:

<a id="eq-nasc-def"></a>

$$
\rho_\text{A} =
\frac{s_\text{A}}{4 \pi \left< \sigma_\text{bs} \right>} =
\frac{s_\text{A}}{\sigma_\text{sp}},
\tag{15}
$$

where $\rho_\text{A}(x)$ is in units of animals nmi<sup>-2</sup>.

[^maclennan_et_al]: MacLennan, D.N., Fernandes, P.G., and Dalen, J. (2002). A consistent approach to definitions and symbols in fisheries acoustics. ICES Journal of Marine Science, 59: 365-369.
[^fisheries_acoustics]: Simmonds, E.J., and MacLennan, D.N. (2005), *Fisheries Acoustics. Theory and Practice.* Blackwell, Oxford. 437 pp.
[^traynor]: Traynor, J.J. (1996). Target-strength measurements of walleye pollock (*Theragra chalcogramma*) and Pacific whiting (*Merluccius productus*). ICES Journal of Marine Science, 53: 253-258.
[^nakken_olsen]: Nakken, O., and Olsen, K. (1977). Target strength measurements of fish. Rapports et Procès-Verbaux des Réunions, 170: 52-69.
[^foote_1986]: Foote, K.G., Aglen, A., and Nakken, O. (1986). Measurement of fish target strength with a split-beam echo sounder. The Journal of the Acoustical Soceity of America, 80: 612-621.
[^foote_1987]: Foote, K.G. (1987). Fish target strengthes for use in echo integrator surveys. Journal of the Acoustical Society of America, 82: 981-987.
[^foote_knudsen_1994]: Foote, K.G., and Knudsen, H.P. (1994). Physical measurement with modern echo integrators. Journal of the Acoustical Society of Japan.