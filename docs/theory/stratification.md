(stratification)=
# Stratification to obtain biomass density


(stratification-intro)=
## Stratification
In practice, the acoustic quantities and biological estimates discussed in previous sections can vary depending on geospatial variation of the biological aggregations themselves. For example, the size and age of fish can vary depending on the survey location, as well as the sex of the fish. Therefore, the acoustic measurements and biological samples are typically stratified to account for these variations, and the biomass density is a function of the stratum (which depends on the geospatial locations) and sex, i.e.

$$
\rho_{B; s} = \rho^i_{B; s}(x,y) = \rho_A(x,y) (\mathbf{L}^i_s)^\top \mathbf{w}^i_s,
$$

where $i$ is the stratum, $\rho_A(x,y)$ is the nautical areal number density at location $(x, y)$, and $\mathbf{L}^i_s$ and $\mathbf{w}^i_s$ are the vectors characterizing the number frequency of fish length and the corresponding weight in stratum $i$, for fish of sex $s$:

$$
\mathbf{L}^i_s = \begin{bmatrix}
L^i_{s,1} \\
L^i_{s,2} \\
L^i_{s,3} \\
\vdots
\end{bmatrix},
$$

and

$$
\mathbf{w}^i_s = \begin{bmatrix}
w^i_{s,1} \\
w^i_{s,2} \\
w^i_{s,3} \\
\vdots
\end{bmatrix}.
$$


Note that the number frequency of fish here is normalized across all length bins and sex, i.e., 

$$
\sum_{s,\ell} L^i_{s,\ell} = 1
$$


### Including age data
In the case when fish age is measured and binned, the biomass density is a function of the stratum (which depends on the geospatial locations), sex, and age:

$$
\rho_{B; s,\alpha} = \rho^i_{B; s,\alpha}(x,y) = \rho_A(x,y) (\mathbf{L}^i_{s,\alpha})^\top \mathbf{w}^i_{s,\alpha},
$$

where $\alpha$ is the age bin,

$$
\mathbf{L}^i_{s,\alpha} = \begin{bmatrix}
L^i_{s,\alpha,1} \\
L^i_{s,\alpha,2} \\
L^i_{s,\alpha,3} \\
\vdots
\end{bmatrix},
$$

and 

$$
\mathbf{w}^i_{s,\alpha} = \begin{bmatrix}
w^i_{s,\alpha,1} \\
w^i_{s,\alpha,2} \\
w^i_{s,\alpha,3} \\
\vdots
\end{bmatrix}.
$$

All of $L^i_{s,\alpha,\ell}$ and $w^i_{s,\alpha,\ell}$ vary depending on the stratum $i$, the fish sex $s$, and the age bin $\alpha$.


Note that the number frequency of fish length here is normalized across all age bins, length bins, and sex within a stratum, i.e.

$$
\sum_{s,\ell,\alpha} L^i_{s,\alpha,\ell} = 1
$$

<!-- ## Jolly-Hampton (1990) stratified sampling 
This analysis provides a coefficient of variation ($\textit{CV}$) for the entire survey by reweighting biomass estimates for $k$ based on the length of each transect $t$ stratified by $i$ to derive estimates of the mean and variance. The first step is summing biomass estimates for each $k$ within each $t$:

$$ 
B^t =
    \sum_{k} B^{k,t}
\label{eq:biomass_transect} \tag{1}
$$

Total transect distance $d_{x,y}^t$ for each $t$ is similarly computed:

$$ 
d_{x,y}^t =
    \sum_{k} d_{x,y}^{k,t}
\label{eq:distance_transect} \tag{2}
$$

These transect distances are then summed for each $i$:

$$ 
D_{x,y} =
    \sum_{t} d_{x,y}^{t}
\label{eq:distance_stratum} \tag{3}
$$

Values of $B^t$ $\eqref{eq:biomass_transect}$, $d_{x,y}^t$ $\eqref{eq:distance_transect}$, and $D_{x,y}$ $\eqref{eq:distance_stratum}$ are then all used to compute the mean transect-length-weighted biomass estimates for each $i$:

$$ 
\tilde{\rho} = 
    \frac{\sum\limits_{t} B^t d_{x,y}^t}{D_{x,y}}
\label{eq:mean_estimate} \tag{4}
$$

Next $B^t$ $\eqref{eq:biomass_transect}$ for each $t$ is standardized to provide a mean biomass-per-distance estimate:

$$ 
\tilde{\rho}^{~t} = 
    \frac{B^t}{d_{x,y}^{t}}
\label{eq:transect_estimate} \tag{5}
$$

which can then be used to compute the squared deviation from the mean along with $\tilde{\rho}$ $\eqref{eq:mean_estimate}$ for each $t$ within $i$:

$$
\tilde{s}^{~t} =
    (\tilde{\rho}^{~t} - \tilde{\rho})^2
\label{eq:squared_deviation} \tag{6}
$$

This can then be used to calculate the sum of weighted squared deviations:

$$
\tilde{s} =
    \sum\limits_{t} w_t^{2} \tilde{s}^t
\label{eq:summed_squared_deviation} \tag{7}
$$

where the stratified weights ($w_t$) for each $t$ are:

$$
w_t = 
    \frac{d_{x,y}^t}{\bar{D}_{x,y}}
\label{eq:stratified_weights} \tag{8}
$$

The variance ($\tilde{\sigma}$) for each $i$ is then calculated:

$$
\tilde{\sigma} =
    \frac{\tilde{s}}{\nu}
\label{eq:variance} \tag{9}
$$

where the $\nu$ represents the degrees of freedom:

$$
\begin{equation}
\nu =
    \begin{cases}
        n^t(n^t-1), & \text{if } n^t > 1 \\
        (n^t)^2, & \text{if } n^t = 1
    \end{cases}
\tag{10} \label{eq:dof} 
\end{equation}
$$

Variance estimates for each $i$ are then weighted by the total area ($A$) for $i$ to compute the overall (weighted) survey variance:

$$
\hat{\sigma} = 
    \sum\tilde{\sigma} A^2
\label{eq:weighted_variance} \tag{11}
$$

where:

$$
A =
    \sum\limits_k A^k
\label{eq:transect_area} \tag{12}
$$

Similar to the weighted variance estimates $\eqref{eq:weighted_variance}$, biomass estimates for each $i$ were also weighted by $A$ to calculate the overall (weighted) survey mean:

$$
\hat{\mu} =
    \sum \tilde{\rho} A
\label{eq:weighted_mean} \tag{13}
$$

The overall survey variance $\eqref{eq:weighted_variance}$ and mean $\eqref{eq:weighted_mean}$ are then both used to calculate $\textit{CV}$:

$$
\textit{CV} =
    \frac{\sqrt{\hat{\sigma}}}{\hat{\mu}}
\label{eq:cv} \tag{14}
$$ -->

## Stratification schemes used in the hake survey
For Pacific hake, two types of stratifications are used:

- **INPFC**: Stratification set by the International North Pacific Fisheries Commission (INFPC) that is based solely on latitude. The US-Canada bienniel hake survey region encompasses 6 strata.
- **KS**: Stratification determined based on the Kolomogorov-Smirnov test for differences of the fish length distributions across survey hauls.
