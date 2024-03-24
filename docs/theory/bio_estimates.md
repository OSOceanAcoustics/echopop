(bio-estimates)=
# Biological estimates


## Number density of scatterers

To obtain the number density of the animal, we define the volume backscattering coefficient ($s_V$, units: m<sup>-1</sup>):

$$
s_V = \rho_V \left< \sigma_{bs} \right>,
$$

and its corresponding logarithmic quantity, the volume backscattering strength ($S_V$, units: dB re 1 m<sup>-1</sup>)

$$
S_V = 10 \log_{10} s_V = 10 \log_{10} \rho_V + 10 \log_{10} \left< \sigma_{bs} \right>
$$

where $\rho_V$ is the number density of scatterers (fish) per unit volume (units: m<sup>-3</sup>).

In fisheries acoustics, we are often interested in quantities per unit area. Therefore, we define the areal backscattering coefficient ($ABC$, or $s_a$, units: m<sup>2</sup>m<sup>-2</sup>)

$$
s_a = s_V H,
$$

where $H$ is the integration height in meter, and the corresponding nautical areal scattering coefficient ($NASC$, or $s_A$, units: m<sup>2</sup>nmi<sup>-2</sup>)

(NASC)=
$$
NASC = s_A = 4 \pi \times 1852^2 \times s_a,
$$

in which the conversion of 1 nmi = 1852 m is used.

Using the above quantities, we obtain

$$
s_a = s_V H = \rho_V \left< \sigma_{bs} \right> H,
$$

Let the areal number density ($\rho_a$, units: m<sup>-2</sup>) be

$$
\rho_a = \rho_V H,
$$

then

$$
s_a = \rho_a \left< \sigma_{bs} \right>.
$$

Similarly, with the corresponding nautical areal number density ($\rho_A$, units: nmi<sup>-2</sup>) being

$$
\rho_A = 1852^2 \rho_a,
$$

then

$$
s_A = NASC = 4 \pi \rho_A \left< \sigma_{bs} \right>.
$$

Note that $NASC$ is the typical output from software packages such as Echoview for biological estimates.





## Biomass estimates

We can obtain an estimate of biomass density ($\rho_B$, units: kg nmi<sup>-2</sup>) by multiplying the areal number density of animals by the average weight ($\left< w \right>$, units: kg)

$$
\rho_B = \rho_A \left< w \right>.
$$

The average weight is

$$
\left< w \right> = \frac{\sum_{j=1}^N w_j}{N},
$$

where $w_j$ is the weight of fish $j$, and $N$ is the total number of fish samples.

In the case when the fish length is binned, which is the case for most fisheries surveys,

$$
\left< w \right> = \mathbf{L}^\top \mathbf{w}.
$$

Here, $\mathbf{L}$ is a vector representing the number frequency $L_\ell$ of fish samples in length bin $\ell$

$$
\mathbf{L} = \begin{bmatrix}
L_1 \\
L_2 \\
L_3 \\
\vdots
\end{bmatrix}
$$

and $\mathbf{w}$ is a vector representing the weight of fish at length $L_\ell$

$$
\mathbf{w} = \begin{bmatrix}
w_1 \\
w_2 \\
w_3 \\
\vdots
\end{bmatrix}.
$$

Note that the number frequency of fish length is normalized across all length bins, i.e.,

$$
\sum_\ell L_\ell = 1.
$$

The $w_\ell$ values can be estimated by the regressed length-weight relationship derived from trawl samples.

With the above quantities, the biomass ($B$, units: kg) can then be estimated by

$$
B = \rho_B A = \rho_A \left< w \right> A,
$$

where $A$ is the unit area associated with the density measure.

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
