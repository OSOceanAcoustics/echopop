(stratification)=
# Stratification to obtain population estimates


(stratification-intro)=
## Stratifying biological length distributions
In practice, the acoustic quantities and biological estimates discussed in previous sections can vary depending on geospatial variation of the biological aggregations themselves. For example, the size and age of fish can vary depending on the survey location, as well as the sex of the fish. Therefore, the acoustic measurements and biological samples are typically stratified to account for these variations, and the biomass density is a function of the stratum (assuming that these are dependent on geospatial locations) and sex:

<a id="intext_eq_29"></a>

(intext_eq_29_md)=
$$
    \rho_{B; s} = \rho^i_{B; s}(x,y) = \rho_A(x,y) (\tilde{\mathbf{L}}^i_s)^\top \mathbf{w}^i_s,
    \tag{2.9}
$$

where $i$ is the stratum, $\rho_A(x,y)$ is the nautical areal number density ({ref}`Eq. 1.15 <eq-115>`) at location $(x, y)$ and $\tilde{\mathbf{L}}^i_s$ is the normalized number frequency ({ref}`Eq. 2.3 <eq-23>` for sex $s$ and stratum $i$):

<a id="intext_eq_210"></a>

(intext_eq_210_md)=
\begin{align*}
    \tilde{\mathbf{L}}^i_s &=
    \left[
        \begin{split}
            \tilde{L}&^i_{s,1} \\ \tilde{L}&^i_{s,2} \\ \tilde{L}&^i_{s,3} \\ &\vdots
        \end{split}
    \right],
    \tag{2.10a}
    \\[2ex]
        \sum_{s,\ell} \tilde{L}^i_{s,\ell} &= 1.
    \tag{2.10b}
\end{align*}

:::{admonition} Defining strata
:class: note
While the descriptions and equations on this page assume that strata ($i$) are geospatial, the same principles apply to other stratification definitions. This can include strata defined using clustering analyses of length distributions or some other statistical unit. It is also assumed that distributions are further separated by sex ($s$).
:::

The weight vector $\mathbf{w}$ ({ref}`Eq. 2.5 <eq-25>`) specific sex $s$ in stratum $i$ ($\mathbf{\tilde{w}}^i_s$) can be normalized similar to $\tilde{\mathbf{L}}$ ({ref}`Eq. 2.3 <eq-23>`; {ref}`Eq. 2.10<eq-210>`). This normalized weight vector is:

<a id="intext_eq_211"></a>

(intext_eq_211_md)=
\begin{align*}
    \tilde{\mathbf{w}}^i_s =
        \left[
        \begin{split}
            \tilde{w}&^i_{s,1} \\ \tilde{w}&^i_{s,2} \\ \tilde{w}&^i_{s,3} \\ &\vdots
        \end{split}
    \right].
    \tag{2.11}
\end{align*}

The quantity $\tilde{w}^i_{s,\ell}$ is defined using the total weights for length $\ell$ ({ref}`Eq. 2.5 <eq-25>`; {ref}`Eq. 2.7 <eq-27>`) specific to sex $s$:

<a id="intext_eq_212"></a>

(intext_eq_212_md)=
\begin{align*}
        \tilde{w}^i_{s,\ell} &= 
            \frac{\mathbf{w}^i_s}{\sum\limits_{\ell} w^i_{s,\ell}} =
                \frac{w^i_\ell}{\sum\limits_{s,\ell} w^i_{s,\ell}},
    \tag{2.12a}    
    \\[2ex]
        \sum\limits_{s,\ell} \tilde{w}^i_{s,\ell} &= 1. 
    \tag{2.12b}
\end{align*}

:::{admonition} Including age ($\alpha$) and other contrasts
:class: note

The above equations can be modified to be binned across any number of contrasts. For instance, adding age ($\alpha$) simplify modifies {ref}`Eq. 2.10 <eq-210>`, {ref}`Eq. 2.11 <eq-211>`, and {ref}`Eq. 2.12 <eq-212>` to:

$$
\begin{align*}
        \tilde{\mathbf{L}}^i_{s,\alpha} = 
            \left[            
                \begin{split}
                    &\tilde{L}^i_{s,\alpha,1} \\ &\tilde{L}^i_{s,\alpha,2} \\ &\tilde{L}^i_{s,\alpha,3} \\ &\hspace{0.5cm} \vdots
                \end{split}
            \right],
\qquad
\sum_{s,\ell,\alpha} \tilde{L}^i_{s,\ell,\alpha} = 1,
\qquad
    \tilde{\mathbf{w}}^i_{s, \alpha} =
        \left[
        \begin{split}
            &\tilde{w}^i_{s,\alpha,1} \\ &\tilde{w}^i_{s,\alpha,2} \\ &\tilde{w}^i_{s,\alpha,3} \\ &\hspace{0.5cm} \vdots
        \end{split}
    \right],
\qquad
\sum_{s,\ell,\alpha} w^i_{s,\ell,\alpha} = 1.
\tag{S.2a--c}
\end{align*}
$$

This would then modify {ref}`Eq. 2.9 <eq-29>` to:

$$
    \rho_{B; s,\alpha} = 
        \rho^i_{B; s,\alpha}(x,y) = \rho_A(x,y) (\tilde{\mathbf{L}}^i_{s,\alpha})^\top \mathbf{\tilde{w}}^i_{s,\alpha},
    \tag{S.2e}
$$

:::

<!-- !!! MOVE TO FEAT IMPLEMENTATION DESCRIPTION -->
<!-- 
## Stratification schemes used in the hake survey
For Pacific hake, two types of stratifications are used:

- **INPFC**: Stratification set by the International North Pacific Fisheries Commission (INFPC) that is based solely on latitude. The US-Canada bienniel hake survey region encompasses 6 strata.
- **KS**: Stratification determined based on the Kolomogorov-Smirnov test for differences of the fish length distributions across survey hauls. -->

(stratified-resampling-algo)=
## Stratified random sampling
{cite:t}`jolly_hampton_1990` proposed a stratified random transect design for estimating the mean and variance of any spatially varying survey quantity. While stratified random sampling as a technique can be generalized to many measured quantities {cite:p}`cochran_sampling_1977`, direct use of the algorithm by {cite:t}`jolly_hampton_1990` requires absolute measurements that are not relative to a unit area (e.g. abundance, biomass). Let $z(\mathbf{x})$ be the surface value at location $\mathbf{x}$ in survey region $D$. In effect, this means that if you are measuring $\rho_\text{A}$, then $z(x^t)$ corresponds is the $\rho_\text{A}$ observed at location $x$ along transect $t$. 

### Estimating the mean

The "true" mean value of these surface values is:

$$
    z = \frac{1}{\|D\|} \int\limits_{D} z(\mathbf{x})\, d\mathbf{x},
    \tag{2.13}
$$

where $\|D\|$ is the total transect length/area.

In practice, each survey transect consists of a sequence of spatial points $\mathbf{x}^k$, where $k$ is the interval (or $\textit{EDSU}$). Here, $k$ is an unique index across the survey. The length of interval $k$, $d^k$, is used to calculate the total distance of transect $t$:

$$
    d^{\,t} = \sum\limits_{k \in t} d^{\,k},
    \tag{2.14}
$$

where the sum is overall intervals $k$ that belong to transect $t$. Quantities like abundance and biomass are not relative to unit areas, so their respective densities can be approximated using $d^{\,t}$:

<a id="intext_eq_215"></a>

(intext_eq_215_md)=
$$
    \hat{z}^{\,t} = 
        \frac{\sum\limits_{k \in t} z(x^k)}{d^t}.
    \tag{2.15}
$$

Per-transect estimates $\hat{z}^{,t}$ are reweighted within each stratum $i$. The weight, $\tau_t$, is based on the ratio between $d^t$ and the average $d^t$ within stratum $i$:

<a id="intext_eq_216"></a>

(intext_eq_216_md)=
$$
    \tau^t = 
        \frac{d^t}{\frac{1}{n^i} \sum\limits_{t \in i} d^t},
    \tag{2.16} 
$$

where $n^i$ is the total number of transects in stratum $i$. The mean density value in stratum $i$ is then:

<a id="intext_eq_217"></a>

(intext_eq_217_md)=
$$
    \hat{z}^{i} = \frac{1}{n^i} \sum\limits_{t \in i} \tau^t \hat{z}^{t}.
    \tag{2.17}
$$

The total survey mean density is then calculated via:

<a id="intext_eq_218"></a>

(intext_eq_218_md)=
$$
    \hat{z} = \frac{\sum\limits_{i} A_i \hat{z}^i}{\sum\limits_i A_i},
    \tag{2.18}
$$

where $A_i$ is the total area (nmi<sup>-2</sup>) of stratum $i$.

### Estimating the variance

The variance of $\hat{z}$, $\mathbb{V}(\hat{z})$, is similarly stratified as the calculation of the mean. The within-stratum variance is defined by:

<a id="intext_eq_219"></a>

(intext_eq_219_md)=
$$
    \mathbb{V}(\hat{z}^i) =
        \frac{\sum\limits_{t \in i} \tau^t (\hat{z}^t - \hat{z}^i)^2}{n^i(n^i - 1)}.
    \tag{2.19}
$$

This is then used to compute the survey variance which is weighted by $A^i$:

<a id="intext_eq_220"></a>

(intext_eq_220_md)=
$$
    \mathbb{V}(\hat{z}) =
        \frac{\sum\limits_i (A^i)^2\, \mathbb{V}(\hat{z}^i)}{\left( \sum\limits_i A^i \right)^2}.
        \tag{2.20}
$$

The coefficient of variation ($\text{CV}$) provides a dimensionless measure of relative uncertainty by relating the total survey variance with its respective mean:

(intext_eq_221_md)=
$$
    \text{CV} = 
        \frac{\sqrt{\mathbb{V}(\hat{z})}}{\hat{z}}.
    \tag{2.21}
$$

(jolly-hampton-bootstrap)=
### Bootstrap resampling and confidence intervals

#### Resampling the estimators

{cite:t}`jolly_hampton_1990` algorithm can also incorporate bootstrapping whereby entire transects are resampled (without replacement) within each stratum, preserving the grouping that defines the spatial process $Z(\mathbf{x})$. For stratum $i$ with $n^i$ transects, each bootstrap replicate $b$ draws a $\hat{p}^*$ proportion of transects:

$$
    m^i_b = 
        \left\lfloor 
            \hat{p}^* n^i + \frac{1}{2}
        \right\rfloor,
    \tag{2.22}    
$$

where $m^i_b$ is the number of resampled transects. Let $\mathbb{T}^*_b$ denote the complese set of resampled transect in bootstrap replicate $b$, where each $t^* \in \mathbb{T}^*_b$ is a transect selected for that replicate. Then the mean density calculation from {ref}`Eq. 2.15 <eq-215>` becomes:

$$
    \hat{z}^{\,t^*}_b = 
        \frac{\sum\limits_{k \, \in \, t^*} z(x^k)}{d^{t^*}}.
    \tag{2.23}
$$

Then the same distance-based weights calculated using {ref}`Eq. 2.16 <eq-216>` is:

\begin{align*}
    \tau^{t^*}_b &= \frac{d^{t^*}}{d^i_b},
    \tag{2.24a}
    \\[2ex]
    d^i_b &= \frac{1}{m^i_b} \sum\limits_{t \, \in \, \mathbb{T}^*_b \, \cap \, i} d^{\,t}.
    \tag{2.24b}
\end{align*}

The stratum mean density calculation from {ref}`Eq. 2.17 <eq-217>` becomes:

$$
    \hat{z}^{i}_b = 
        \frac{1}{m^i_b} \sum\limits_{t \, \in \, \mathbb{T}^*_b \, \cap \, i} \tau^{t} \hat{z}^{t}_b.
    \tag{2.25}
$$


The stratum replicates are then pooled across strata via the area-weighted mean in {ref}`Eq. 2.18 <eq-218>` to obtain a survey-level replicate $\hat{z}_b$:

$$
    \hat{z}_b = \frac{\sum\limits_{i} A^i \hat{z}^i_b}{\sum\limits_i A^i},
    \tag{2.16}
$$

where $A^i_b$ is the coverage area corresponding to the resampled transects $t^* \in \mathbb{T}^*_b$ for replicate $b$.

Similarly, the resampled variance for stratum $i$ ({ref}`Eq. 2.19 <eq-219>`) and survey region $D$ ({ref}`Eq. 2.20 <eq-220>`) per replicable $b$ is expressed as:


\begin{align*}
    \mathbb{V}(\hat{z}^i_b) &=
        \frac{\sum\limits_{t \, \in \, \mathbb{T}^*_b \, \cap \, i} \tau^t_b (\hat{z}^t_b - \hat{z}^i_b)^2}{m^i(m^i - 1)},
    \tag{2.27a}
    \\[2ex]
    \mathbb{V}(\hat{z}_b) &= 
        \frac{\sum\limits_i (A^i_b)^2\, \mathbb{V}(\hat{z}^i)}{\left( \sum\limits_i A^i_b \right)^2}.
    \tag{2.27b}
\end{align*}

Finally, the mean density and variance for each $b$ are used to compute the replicate $\text{CV}$:

$$
    \text{CV}_b = 
        \frac{\sqrt{\mathbb{V}(\hat{z}_b)}}{\hat{z}_b}.
    \tag{2.27}
$$

(jolly-hampton-ci)=
#### Confidence interval estimation

After generating bootstrap replicates for the $\hat{z}_b$, $\mathbb{V}(\hat{z}_b)$variance, and $\text{CV}$, confidence intervals for these quantities can be constructed using a variety of methods. One common approach is using the **percentile method** {cite:p}`efron_bootstrap_1994` that has historically been widely used in survey statistics {cite:p}`cochran_sampling_1977`. Let $\{\hat{z}_b\}_{b=1}^B$ denote the set of bootstrap replicates for the survey mean density, where $B$ is the total number of bootstrap samples. The confidence interval for the survey mean at level $(1-\alpha)$ is:

$$
    \left[
        \hat{z}_{(\alpha/2)},\;
        \hat{z}_{(1-\alpha/2)}
    \right],
    \tag{2.28}
$$

where $\hat{z}_{(q)}$ is the $q$-th quantile of the sorted bootstrap replicates. This method can be applied analogously to other survey statistics, such as the $\text{CV}_b$ or stratum-level means and variances.

The **empirical distribution** (also commonly referred to as the **delta method**) is another method that can account for resampling bias and is particularly useful when the bootstrap replicates not symmetrically distributed around the population statistic {cite:p}`efron_nonparametric_1981`. Given bootstrap replicates $\{z_b\}_{b=1}^B$ and the original population statistic $\hat{\theta}$, the deviation for each replicate is:

$$
    \delta_b = \hat{z}_b - \hat{\theta}.
    \tag{2.29}
$$

For a confidence level $(1-\alpha)$, let $q_1 = \alpha/2$ and $q_2 = 1-\alpha/2$. The empirical confidence interval is then:

$$
    \left[
        \bar{\hat{z}_b} + Q_{q_1}(\delta),\;
        \bar{\hat{z}_b} + Q_{q_2}(\delta)
    \right],
    \tag{2.30}
$$

where $\bar{\hat{z}_b}$ is the mean of the bootstrap replicates and $Q_{q}(\delta)$ is the $q$-th quantile of the deviation distribution.