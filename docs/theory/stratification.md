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






## Stratification schemes used in the hake survey
For Pacific hake, two types of stratifications are used:

- **INPFC**: Stratification set by the International North Pacific Fisheries Commission (INFPC) that is based solely on latitude. The US-Canada bienniel hake survey region encompasses 6 strata.
- **KS**: Stratification determined based on the Kolomogorov-Smirnov test for differences of the fish length distributions across survey hauls.






## Jolly-Hampton stratified sampling 
Jolly and Hampton {cite:p}`jolly_hampton_1990` proposed a method to provide a coefficient of variation ($\textit{CV}$) for the entire survey by weighting biomass estimates across all strata to derive estimates of the mean and variance of the biomass density. 

Given the biomass $B^k$ for all intervals $k$ in transect $t$ within stratum $i$, the biomass density per unit area is

$$
\hat{\rho}_B^t = \frac{ B^t }{ d^t },
$$
where $B^t=\sum_{k \in t} B^k$ is the total biomass and $d^t=\sum_{k \in t} d^k$ is the total distance of transect $t$.

The mean weighted biomass density in stratum $i$ is then

$$ 
\hat{\rho}_B^i = \sum_{t \in i} \gamma_t B^t,
$$

where $\gamma_t = d^t / d^i$ is the weight of transect $t$, and $d^i = \sum_{t \in i} d^t$ is the total transect distance in stratum $i$.

Across the entire survey area, the biomass density estimate can similarly be obtained by a weighted mean over all strata:

$$
\hat{\rho}_B = \frac{ \sum_i A_i \hat{\rho}_B^i }{ \sum_i A_i },
$$

where $A_i$ is the area of stratum $i$.

 
When transects are randomly selected within the strata, the variance of the total biomass can be estimated as

$$
\textrm{Var}(\hat{\rho}_B) = \frac{ \sum_i A_i^2 \textrm{Var}(\hat{\rho}_B^i) }{ \left( \sum_i A_i \right)^2 },
$$

where

$$
\textrm{Var}(\hat{\rho}_B^i) = \frac{ \sum_{t \in i} \gamma_t^2 \left( \hat{\rho}_B^t - \hat{\rho}_B^i \right)^2 }{ n^i (n^i-1) },
$$

in which $n^i$ is the number of transects in stratum $i$.

The $\textrm{CV}$ can then be calculated using:

$$
CV = \frac{ \sqrt{\textrm{Var}(\hat{\rho}_B)} }{\hat{\rho}_B}.
$$
