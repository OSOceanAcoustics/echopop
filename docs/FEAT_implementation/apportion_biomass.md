(apportion-biomass)=
# Apportioning kriged biomass density

(apportion-aged)=
## Apportioning within aged fish samples

The weight proportions of aged fish samples of sex $s$, length $\ell$, and age $\alpha$ with respect to all aged _and_ unaged fish samples are used to apportion the aged component of the total biomass. The total weight for stratum $i$ is:

$$
    w_i^{\text{all}} = w_i^{\text{aged}} + w_i^\text{unaged},
$$

which can then calculate the within-group weight proportions for aged fish:

$$
    {r_w}^\textrm{aged/all}_{i, s,\ell,\alpha} = 
    \frac{ \sum_{j \in J^\textrm{aged}_{i,s,\ell,\alpha}} w_j }{ w_i^{\text{all}}  }.
$$

This proportion can then be applied to the transect-based estimates to apportion biomass for stratum, length, age, and sex. The apportioned biomass for transect interval $k$ is therefore:

$$
    B^{k, \textrm{aged}}_{i, s,\ell,\alpha} = 
    \rho^k_\text{B} A^k {r_w}^\textrm{aged/all}_{i, s,\ell,\alpha},
$$

where $A^k$ is the area of the $k$<sup>th</sup> interval.

The apportioned biomass over all transect intervals for sex $s$ is then:

$$
    B^\textrm{aged}_{s,\ell,\alpha} =
        \sum_i \sum_k B^{k, \textrm{aged}}_{i, s,\ell,\alpha}.
$$

It follows naturally that the combined biomass quantities for male and female fish are:

$$
B^{k, \textrm{aged}}_{i, \ell,\alpha} = B^{k, \textrm{aged}}_{i,\textrm{M},\ell,\alpha} + B^{k, \textrm{aged}}_{i,\textrm{F},\ell,\alpha}.
$$

(apportion-unaged)=
## Apportioning within unaged fish sampless

For unaged fish, the procedure is more complicated, because the samples are not weighed individually, but are weighed only for each sex separately. Therefore, the length-weight relationship ($\mathcal{W}(\ell)$) derived from aged fish samples ({ref}`Eq. (2.6) <eq-26>`, {ref}`Eq. (2.7) <eq-27>`) is used to partition the total weight into different length bins. Within each length bin, the age distribution of aged fish samples of that length is further used to partition the weight into different length-age bins.

### Unaged weight proportions of different length bins

The number proportion of unaged fish samples of length $\ell$ and sex $s$ with respect to all unaged fish samples is:

$$
    {r_n}^\textrm{unaged/unaged}_{i,s,\ell} = \frac{ n^\textrm{unaged}_{i,s,\ell} }{ n^\textrm{unaged}_i },
$$

where $n^\textrm{unaged}_i$ is the total number of fish in stratum $i$. The proportion can also be integrated over $s$ via:

$$
    {r_n}^\textrm{unaged/unaged}_{i,\ell} = \sum_s {r_n}^\textrm{unaged/unaged}_{i,s,\ell}.
$$

It is worth noting that these quantities are analagous to $\tilde{\mathbf{L}}$ described by {ref}`Eq. (2.3) <eq-23>`. Using the $\mathcal{W}(\ell)$ computed for all aged fish samples (i.e., inclusive of both male and female fish), the weight proportion of unaged fish samples of length $\ell$ with respect to all unaged fish samples is:

$$
    {r_w}^\textrm{unaged/unaged}_{i,\ell} =
        \frac{
                {r_n}^\textrm{unaged/unaged}_{i,\ell} \times \mathcal{W}(\ell)_\textrm{all}
            }{
                \sum\limits_{\ell} {r_n}^\textrm{unaged/unaged}_{i,\ell} \times \mathcal{W}(\ell)_\textrm{all}
            } .
$$

### Unaged weight proportions of different sexes
Separately, the weight proportions of male and female within the unaged fish samples can similarly be calculated by using the number proportion of a given sex and the sex-specific length-weight relationship.

The inferred total weights of male and female unaged fish are:

$$
\begin{align*}
    w^\textrm{unaged}_\textrm{i,M} &= \sum_\ell {r_n}^\textrm{unaged/unaged}_{i,\textrm{M},\ell} \, n^\textrm{unaged}_\textrm{i,M} \, \mathcal{W}(\ell)_\textrm{all} \\
    w^\textrm{unaged}_\textrm{i,F} &= \sum_\ell {r_n}^\textrm{unaged/unaged}_{i, \textrm{F},\ell} \, n^\textrm{unaged}_\textrm{i,F} \, \mathcal{W}(\ell)_\textrm{all} .
\end{align*}
$$

The weight proportions of male and female unaged fish with respect to all unaged fish samples are then:

$$
\begin{align*}
    {r_w}^\textrm{unaged/unaged}_\textrm{i,M} &= 
    \frac{ w^\textrm{unaged}_\textrm{i,M} }{ w^\textrm{unaged}_\textrm{M} + w^\textrm{unaged}_\textrm{F} } \\
    {r_w}^\textrm{unaged/unaged}_\textrm{i,F} &= \frac{ w^\textrm{unaged}_\textrm{i,M} }{ w^\textrm{unaged}_\textrm{i,M} + w^\textrm{unaged}_\textrm{i,F} } .
\end{align*}
$$

The weight proportions of unaged male and female fish with respect to all fish samples (aged _and_ unaged) are then:

$$
\begin{align*}
    {r_w}^\textrm{unaged/all}_\textrm{i,M} &= {r_w}^\textrm{unaged/unaged}_\textrm{i,M} \, {r_w}^\textrm{unaged/all}_i \\
    {r_w}^\textrm{unaged/all}_\textrm{i,F} &= {r_w}^\textrm{unaged/unaged}_\textrm{i,M} \, {r_w}^\textrm{unaged/all}_i .
\end{align*}
$$

### Unaged biomass apportioned with sex and length

The sex-specific and overall biomass of unaged fish in interval $k$ can then be expressed as:

$$
\begin{align*}
    B^{k, \textrm{unaged}}_{i, s, \ell} &= \rho^k_\text{B} \, A^k \, {r_w}^\textrm{unaged/all}_{i,s} \, {r_w}^\textrm{unaged/unaged}_{i,\ell} \\
    B^{k, \textrm{unaged}}_{i, \ell} &= \rho^k_\text{B} \, A^k \, {r_w}^\textrm{unaged/all}_{i} \, {r_w}^\textrm{unaged/unaged}_{i,\ell}.
\end{align*}
$$

Summing across all transect intervals, the sex-specific and overall apportioned biomass of unaged fish is then:

$$
\begin{align*}
    B^\textrm{unaged}_{i, s, \ell} &= \sum_k B^{k, \textrm{unaged}}_{i, s, \ell} \\
    B^\textrm{unaged}_{i, \ell} &= \sum_k B^{k, \textrm{unaged}}_{i, \ell}.
\end{align*}
$$

(unaged-biomass-apportionment)=
### Unaged biomass apportioned with sex, length, and age

To further partition the biomass of different length bins into different age bins, the weight proportion of fish samples of age $\alpha$ with respect to all aged fish samples of length $\ell$ is used. The calculation is carried out for male and female fish separately before the estimates are combined. The total apportioned biomass of unaged male and female fish of length $\ell$ and age $\alpha$ is:

$$
\begin{align*}
    B^\textrm{unaged}_{i, \textrm{M},\ell,\alpha} &= B^\textrm{unaged}_{i,\textrm{M},\ell} \times \frac{ B^\textrm{aged}_{i,\textrm{M},\ell,\alpha} }{ \sum_\alpha B^\textrm{aged}_{i,\textrm{M},\ell,\alpha} } \\
    B^\textrm{unaged}_{i,\textrm{F},\ell,\alpha} &= B^\textrm{unaged}_{i,\textrm{F},\ell} \times \frac{ B^\textrm{aged}_{i,\textrm{F},\ell,\alpha} }{ \sum_\alpha B^\textrm{aged}_{i,\textrm{F},\ell,\alpha} } ,
\end{align*}
$$

where the total apportioned biomass of unaged fish across length $\ell$ and age $\alpha$ is expressed as:

$$
    B^\textrm{unaged}_{i,\ell,\alpha} = B^\textrm{unaged}_{i,\textrm{M},\ell,\alpha} + B^\textrm{unaged}_{i,\textrm{F},\ell,\alpha}.
$$

## Total biomass apportioned with length and age

Using the quantities from above, the total apportioned biomass of all (aged _and_ unaged) fish of length $\ell$ and age $\alpha$ is then

$$
    B^\textrm{all}_{i,\ell,\alpha} = B^\textrm{aged}_{i,\ell,\alpha} + B^\textrm{unaged}_{i,\ell,\alpha}.
$$
