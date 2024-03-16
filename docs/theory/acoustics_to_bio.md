# Acoustic backscatter to biological estimates


## TS-length regression

$$
\text{TS} = mL + b
$$

where $\text{TS}$ is the target strength (units: dB re. 1 m<sup>-2</sup>), $L$ is the total length (units: m), $m$ is the slope, and $b$ is the <i>y</i>-intercept.

Here, the backscattering cross-section of a given animal $\sigma_{bs}$ is related to $\text{TS}$ by

$$
\text{TS} = 10 \log_{10} \sigma_\mathrm{bs}.
$$

For a group of animal in stratum $i$ with different length (and hence different $\text{TS}$ and $\sigma_\mathrm{bs}$), the mean backscattering cross-section is

$$
\bar{\sigma}_\mathrm{bs}^i = \frac{\sum_{j=1}^n \sigma_{\mathrm{bs}, j}^i}{n}.
$$

<!-- ## Imputation

Let $\hat{i}$ represents the expected strata, $\hat{i}_{\mathrm{miss}} = i$, and $ \hat{i}$ which represents values of $i$ missing from $\hat{i}$

$$
\bar{\sigma}_{\mathrm{bs}}^{i} = \begin{cases}
    \bar{\sigma}_{\mathrm{bs}}^{i+1} & \text{if } i = \hat{i}_{\mathrm{min}}  \text{ and } i + 1 \in \hat{i} \\
    \bar{\sigma}_{\mathrm{bs}}^{i-1} & \text{if } i = \hat{i}_{\mathrm{max}}  \text{ and } i + 1 \in \hat{i} \\
    \frac{1}{2}(\bar{\sigma}_{\mathrm{bs}}^{i-1} + \bar{\sigma}_{\mathrm{bs}}^{i+1}) & \text{if } i \in \hat{i}_{\mathrm{miss}} \text{ and } (i-1, i+1) \subseteq \hat{i} \\
    \bar{\sigma}_{\mathrm{bs}}^{i} & \text{if } i \in \hat{i} 
\end{cases}
$$ -->


## Converting NASC to biological estimates

The areal numerical density $\rho_N^k$ (units: nmi<sup>-2</sup>) of a target animal species with an averaged backscattering cross section $\bar{\sigma}^i_{bs}$ can be obtained by:

$$
\rho_N^k = \frac{\text{NASC}^k}{4 \pi \bar{\sigma}_{\mathrm{bs}}^i} ,
$$

where $\text{NASC}^k$ is the nautical areal scattering coefficient in the along-transect segment $k$, with the transect being in stratum $i$.

To convert the areal number density to areal _biomass_ density $\rho_W^k$: 

$$
\rho_W^k = \rho_N^k W^i ,
$$

where $W^i$ is the length-weight conversion in stratum $i$.

Using the above quantities, for the along-transect segment $k$, the animal abundance $N$:

$$
N^k = \rho_N^k A^k . 
$$

The biomass $W^k$ is:

$$
W^k = \rho_W^k A^k
$$

And the apportioned biomass is:

$$
W_{\alpha, s, \ell}^k = W^k p_{\alpha, s, \ell}^i ,
$$

where $A^k$ is the area of along-transect segment $k$, and $p^i_{\alpha, s, \ell}$ is the fraction of the animal of age $\alpha$, length $\ell$, and sex $s$ in stratum $i$. Recall that the along-transect segment $k$ is in stratum $i$.