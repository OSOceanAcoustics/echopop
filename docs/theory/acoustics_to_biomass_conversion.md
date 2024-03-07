# Acoustic backscatter to biomass conversion

![ text ](../images/example_indexing.jpg)

TS-length regression

$$\textit{TS} = mL_{\textit{cm}} + b$$

where $\textit{TS}$ is the target strength (dB re. 1 m<sup>-2</sup>), $L_{\textit{cm}}$ is the total length (cm), $m$ is the slope, and $b$ is the <i>y</i>-intercept.

$$\hat{\sigma_\mathrm{bs}} = 10^{\textit{TS}(L)/10}$$

$$\bar{\sigma}_{\mathrm{bs}}^{i} = \frac{\sum\limits_{i=0} \hat{\sigma}_{\mathrm{bs}}^{i}(L)}{n^{i}}$$

Imputation

$\hat{i}$ represents the expected strata. 
$\hat{i}_{\mathrm{miss}} = i \ \hat{i}$ which represents values of $i$ missing from $\hat{i}$

$$
\bar{\sigma}_{\mathrm{bs}}^{i} = \begin{cases}
    \bar{\sigma}_{\mathrm{bs}}^{i+1} & \text{if } i = \hat{i}_{\mathrm{min}}  \text{ and } i + 1 \in \hat{i} \\
    \bar{\sigma}_{\mathrm{bs}}^{i-1} & \text{if } i = \hat{i}_{\mathrm{max}}  \text{ and } i + 1 \in \hat{i} \\
    \frac{1}{2}(\bar{\sigma}_{\mathrm{bs}}^{i-1} + \bar{\sigma}_{\mathrm{bs}}^{i+1}) & \text{if } i \in \hat{i}_{\mathrm{miss}} \text{ and } (i-1, i+1) \subseteq \hat{i} \\
    \bar{\sigma}_{\mathrm{bs}}^{i} & \text{if } i \in \hat{i} 
\end{cases}
$$
NASC to areal number density:

$$\hat{\rho}_{A,N}^{i,j,k} = \frac{\textit{NASC}^{i,j,k}}{4 \pi \bar{\sigma}_{\mathrm{bs}}^{i}}$$

where $\hat{\rho}_{A,N}$ is the areal number density (animals nmi<sup>-2</sup>) or abundance of animals ($N$) per unit area ($A$). The subscript indices represent the 
stratum ($i$), transect ($j$), and along-transect segment ($k$).

Areal number to biomass density:

$$\hat{\rho}_{A,B}^{i,j,k} = \hat{\rho}_{A,N}^{i,j,k} \hat{w}^{i}$$

where $\hat{w}^{i}$ is the derived length-weight conversion in stratum $i$.

Abundance: 

$$N^{i,j,k} = \hat{\rho}_{A,N}^{i,j,k} A^{i,j,k}$$

Biomass: 

$$B^{i,j,k} = \hat{\rho}_{A,B}^{i,j,k} A^{i,j,k}$$

Indexed biomass:

$$B_{\alpha, s, \ell}^{i,j,k} = B^{i,j,k} p_{\alpha, s, \ell}^{i,j,k}$$

where $\alpha$ is binned age, $s$ is sex, $\ell$ is the binned length, and $p$ is the relative proportion.