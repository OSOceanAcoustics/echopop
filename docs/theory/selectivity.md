(net-selectivity)=
# Correcting for net-selectivity

This section details the theoretical and mathematical framework used to adjust observed catch proportions to estimate the true population structure.

---

## Gear selectivity
Fishery surveys are subject to **gear selectivity**, where the probability of a fish being retained by the gear is a function of its physical dimensions, which is primarily body length. Trawl nets act as a "length-based filter" where smaller fish may pass through the mesh while larger fish are retained. Because the sampling gear is biased by size, the resulting length- and length-at-age distributions are also biased. The inverse of the gear's selectivity estimates the population scale by re-weighting the observed distributions. This reconstruction follows the standard methodology for measuring and correcting the selectivity of towed fishing gears [^1].

---

## The logistic selectivity model

The logistic function, also known as a **selection ogive**, is a standard mathematical model for describing the retention of fish by towed sampling gears. It belongs to the class of **sigmoidal functions**, characterized by a distinct "S-shaped" profile that maps a real-valued input (length) to a probability space between 0 and 1. The logistic model is favored because it mirrors the cumulative probability of two distinct processes:

1. **Physical encounter**: the probability a fish of length $L$ physically contacts the gear.
2. **Retention event**: the probability that, once inside the net, the fish is unable to pass through the mesh.

The **sigmoidal** shape of the curve effectively captures the transition from **complete escapement** (very small fish) to **complete retention** (very large fish) [^3]. The symmetry of the logistic curve assumes that the probability of escape decreases at the same rate that the probability of retention increases as fish size approaches the $L_{50}$ threshold. This threshold is the length at which a fish has a 50% chance of being retained.

### Mathematical derivation of coefficients

When using a logistic selectivity model, parameters like $L_{50}$ are derived from the the slope ($\beta_1$) and intercept ($\beta_0$) of the regression. By applying a logit transformation to the retention probability, $S(L)$, we move from a bounded S-curve to a linear equation:

$$
    \ln\left( \frac{S(L)}{1 - S(L)} \right) = \beta_0 + \beta_1 L,
    \tag{2.31}
$$

where $\beta_0$ represents the log-dds of retention when length is zero and $\beta_1$ is the change in log-odds of retention for every unit increase in length.

### Biological interpretability

To make these coefficients biologically interpretable, they are transformed into the length at 50% retention ($L_{50}$) and the selection range ($SR$). The $L_{50}$ occurs where the probability of retention is exactly $0.5$ (where the log-odds are zero):

$$
    0 = \beta_0 + \beta_1 L_{50} \implies L_{50} = -\frac{\beta_0}{\beta_1}.
    \tag{2.32}
$$

Similarly, the selection range ($SR$) is defined as the length interval between $25\%$ and $75\%$ retention ($L_{75} - L_{25}$), which is inversely proportion to the steepness of the slope $\beta_1$:

$$
    SR = \frac{2 \ln(3)}{\beta_1}.
    \tag{2.33}
$$

### Parameterization

The selectivity probability $S(L)$ can be expressed using either the regression-based ($\beta_0$, $\beta_1$) or biological parameters ($L_{50}$, $SR$):

\begin{align*}
        S(L) &= \frac{1}{1 + e^{-(\beta_0 + \beta_1 L)}},
    \tag{2.34a}
    \\[2ex]
        S(L) &= \left[ 1 + e^{\left( \frac{2 \ln(3) (L_{50} - L)}{SR} \right)} \right]^{-1}.
    \tag{2.34b}
\end{align*}

---

## Correcting number proportions
For a [distribution consisting of only length](./bio_estimates.md#intext_eq_24_md) ($\tilde{L}_\ell$), the adjusted population proportion $\hat{L}_\ell$ is calculated by dividing the observed proportion $\tilde{L}_\ell$ by its corresponding selectivity $S(L_\ell)$ [^2]:

\begin{align*}
        \hat{L}_\ell &= \frac{\tilde{L}_\ell S(L_\ell)^{-1}}{\sum\limits_{i} (\tilde{L}_i S(L_i)^{-1})},
    \tag{2.35a}
    \\[2ex]
        \sum\limits_\ell \hat{L}_\ell &= 1,
    \tag{2.35b}
\end{align*}

where $\hat{L}_\ell$ is re-normalized such that it sums to 1.

For a [distribution consisting of both length and age](./stratification.md#indexing-by-length-age) ($\tilde{L}_{\ell, \alpha}$), the selectivity remains strictly a function of length. The 1D selectivity vector is broadcast across the age dimension, applying the same retention probability to all age classes within a specific length bin:

\begin{align*}
        \hat{L}_{\ell, \alpha} &= \frac{\tilde{L}_{\ell, \alpha} S(L_{\ell, \alpha})^{-1}}{\sum\limits_{i, j} (\tilde{L}_{i,j} S(L_{i,j})^{-1})},
    \tag{2.36a}
    \\[2ex]
        \sum\limits_{\ell , \alpha}\hat{L}_{\ell, \alpha} &= 1.
    \tag{2.36b}
\end{align*}

The length-based selectivity is applied uniformly over age due to the assumption that a net's physical mesh interacts with the morphology of a fish (i.e., length and girth). Within a single bin $\ell$, individuals are assumed to have the same probability of capture regardless of their age $a$ or sex $s$. However, because younger fish are statistically more likely to occupy smaller length bins with lower $S(L)$, this length-based weighting naturally shifts the resulting marginal age distribution toward younger cohorts, correcting the "missing" data.

---

## Practical considerations and constraints

### Minimum Selectivity ($S_\text{min}$)
In the extreme left tail of the selection ogive, $S(L)$ approaches zero. Dividing by near-zero values can lead to extreme variance inflation. To stabilize the reconstruction, a lower bound is enforced:

$$
    S_\text{eff}(L) = \max(S(L), S_\text{min}).
    \tag{2.37}
$$

### Zero-Catch Handling
The correction is a linear weighting operator. If the observed proportion for a bin is $0$, the adjusted proportion remains $0$. This assumes that a zero catch in a specific stratum represents a true absence or a sampling artifact that cannot be reconstructed without external prior information. Additional methods are required to account for potentially "false" zeros in the observed distributions.

---

## Footnotes

[^1]: Wileman, D. A., Ferro, R. S. T., Fonteyne, R., & Millar, R. B. (Eds.). (1996). *Manual of methods of measuring the selectivity of towed fishing gears*. ICES Cooperative Research Report No. 215.
[^2]: Millar, R. B., & Fryer, R. J. (1999). *Estimating the size-selection curves of towed selection, snapshot, and multi-mesh pelagic sampling gears*. Reviews in Fish Biology and Fisheries, 9(1), 89-116.
[^3]: Sparre, P., & Venema, S. C. (1998). *Introduction to tropical fish stock assessment. Part 1: Manual*. FAO Fisheries Technical Paper No. 306.1, Rev. 2.
