(acoustics-basics)=
# Acoustics basics

## Backscattering cross-section and target strength
For a given scatterer, the differential backscattering cross section ($\sigma_{bs}$, units: m<sup>2</sup>) and backscattering cross section ($\sigma_{b}$, units: m<sup>2</sup>) are related by

$$
\sigma_{bs} = \frac{\sigma_b}{ 4 \pi}
$$

The target strength ($TS$, units: dB re 1 m<sup>2</sup>) is defined as

$$
TS = 10 \log_{10} \sigma_{bs}
$$

For a group of $N$ animals, the mean differential backscattering cross-section is

$$
\left< \sigma_{bs} \right> = \frac{\sum_{j=1}^N \sigma_{bs,j} }{ N },
$$

where $\sigma_{bs,j}$ is the differential backscattering cross-section of animal $j$, which often varies as a function of its length $L_j$:

$$
\sigma_{bs,j} = \sigma_{bs,j}(L_j)
$$




## TS-length relationship

One common avenute to estimate TS of a scatterer is based on empirical relationshpi between TS and length, which can be expressed by

$$
TS = mL + b,
$$

where $L$ is the total length, $m$ is the slope, and $b$ is the <i>y</i>-intercept.

For Pacific hake, the empitical TS-length relationship used in Echopop is

$$
TS = 20 L - 68
$$

where $L$ is the fish fork length in cm.

Therefore, for Pacific hake

$$
\sigma_{bs} = 10^{TS/10} = 10^{-6.8} L^2
$$
