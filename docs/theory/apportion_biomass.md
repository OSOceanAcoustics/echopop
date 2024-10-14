(apportion-biomass)=
# Apportioning kriged biomass density

The challenges associated with apportioning the kriged biomass is to properly combine the various pieces of information to distribute the biomass into different sex, length, and age groups. This is because fish samples obtained from a haul (trawl) are processed at two different stations that report different biometric data:

(apportion-trawl-stations)=
- Station 1 (the "unaged" fish samples):
    - The fish samples are first sorted according to sex. The length of fish samples of the same sex are measured and binned into a length distribution. Only the total weight of all unaged fish samples from a haul is reported.
    - `length_df`: contains the length distribution of unaged fish of different sexes
    - `haul_df`: contains the total weight of all unaged fish samples in each haul
- Station 2 (the "aged" fish samples):
    - The sex, length, and weight of each fish sample are measured at sea, and the age of each fish samples is measured in the lab on shore after the survey.
    - `specimen_df`: contains the sex, length, weight, and age information for each aged fish sample

Therefore, to apportion the total biomass into different sexes, length bins, and age bins, weight of the aged fish samples can be used directly (see [](apportion-aged)), but weight of the unaged fish samples needs to be further partitioned according to information obtained from the aged fish samples (see [](apportion-unaged)).





## Symbols and notation

### Transect-related parameters

- $i$ denotes a specific stratum
- $h$ denotes a specific haul, which belongs to a particular stratum
- $k$ denotes a specific transect interval, which belongs to a particular stratum

These transec-related parameters will appear as superscript of a given quantity.

### Biometric parameters
- $s$ is the sex of fish sample and can take values $\textrm{M}$ (male) or $\textrm{F}$ (female)
- $\ell$ is the length of fish sample
- $\alpha$ is the age of fish sample

These biometric parameters will appear as subscript of a given quantity.

### Aged and unaged samples
A superscript is used to denote the set of fish samples:
- $\textrm{all}$: all aged and unaged fish samples
- $\textrm{aged}$: all aged fish samples
- $\textrm{unaged}$: all unaged fish samples

### Quantities
- $w_j$ is the weight of fish sample $j$
- $J$ denotes a set of fish samples. For example:
    - $J^\textrm{all}$ is the set of all aged _and_ unaged fish samples
    - $J^{\textrm{unaged}, h}$ is the set of unaged fish samples from haul $h$
    - $J^\textrm{aged}_{s,\ell,\alpha}$ is the set of aged fish samples of sex $s$, length $\ell$, and age $\alpha$
- $n$ denotes the number of fish samples. For example:
    - $n^\textrm{unaged}$ is the total number of unaged fish samples
    - $n^\textrm{unaged}_\ell$ is the number of unaged fish samples of length $\ell$
- $W_s(\ell)$ describes the length-weight relationship for a population of fish of sex $s$. For example:
    - $W_\textrm{M}(\ell)$ is the length-weight relationship for male fish
    - $W(\ell)$ is the length-weight relationship for all fish (male and female combined)
- $r_w$ denotes the weight proportion
    - ${r_w}^\textrm{aged/all}$ is the weight proportion of aged fish with reference to all fish samples
    - ${r_w}^\textrm{aged/aged}_{s,\ell}$ is the weight propotion of aged fish of sex $s$ and length $\ell$ with reference to all aged fish samples
- $\rho^k_B$ is the biomass density in transect interval $k$
- $A^k$ is the area covered by transect interval $k$
- $B^k$ is the biomass in transect interval $k$

```{attention}
The sex, length, and age distributions of fish vary across strata (see [](stratification-intro)). Therefore, the specific strata is implicit in the description below. Most quantities can be annotated with an additional superscript $i$ to indicate the stratum the sample or quantity belongs to. Exceptions arise when the quantity is specific for a transect interval ($k$) or a haul ($k$), since a transect interval or a haul only belongs to one stratum.
```




## Overall weight proportion of aged and unaged fish samples

### Aged fish samples
For aged fish samples, the data is given in the form of individual fish samples with sex, length, and age information (in `specimen_df`). Therefore, the summed weight is computed by summing over both male and female aged fish samples of all length and age bins:

$$
\begin{split}
w^\textrm{aged}
&= w^\textrm{aged}_\textrm{M} + w^\textrm{aged}_\textrm{F} \\
&= \sum_{\ell,\alpha} w^\textrm{aged}_{\textrm{M},\ell,\alpha} + \sum_{\ell,\alpha} w^\textrm{aged}_{\textrm{F},\ell,\alpha} \\
&= \sum_{\ell,\alpha} \sum_{j \in J^\textrm{aged}_{\textrm{M},\ell,\alpha}} w_j + \sum_{\ell,\alpha} \sum_{j \in J^\textrm{aged}_{\textrm{F},\ell,\alpha}} w_j
\end{split}
$$

### Unaged fish samples
For unaged fish samples, the data is given in the form of haul catch weight (in `haul_df`) and no individual weight information is available. Therefore, the total summed weight of unaged fish is

$$
w^\textrm{unaged} = \sum_{j \in J^{\textrm{unaged},h}} w_j
$$

### Overall weight proportion
The overall weight propotion of aged and unaged fish are then

$$
{r_w}^\textrm{aged/all} = \frac{ w^\textrm{aged} }{ w^\textrm{aged} + w^\textrm{unaged} } 
$$

and

$$
{r_w}^\textrm{unaged/all} = \frac{ w^\textrm{unaged} }{ w^\textrm{aged} + w^\textrm{unaged} } 
$$

respectively.





(apportion-aged)=
## Apportioning within aged fish samples

The weight proportions of aged fish samples of sex $s$, length $\ell$, and age $\alpha$ with respect to all aged _and_ unaged fish samples are used to apportion the aged component of the total biomass.

The weight propotions are:

$$
{r_w}^\textrm{aged/all}_{s,\ell,\alpha} = \frac{ \sum_{j \in J^\textrm{aged}_{s,\ell,\alpha}} w_j }{ \sum_{j \in J^\textrm{all}} w_j }
$$

The apportioned biomass for transect interval $k$ is

$$
B^{k, \textrm{aged}}_{s,\ell,\alpha} = \rho^k_B A^k {r_w}^\textrm{aged/all}_{s,\ell,\alpha}
$$

The apportioned biomass over all transect intervals is then

$$
B^\textrm{aged}_{s,\ell,\alpha} = \sum_k B^{k, \textrm{aged}}_{s,\ell,\alpha}
$$

It follows naturally that the combined biomass quantities for male and female fish are:

$$
B^{k, \textrm{aged}}_{\ell,\alpha} = B^{k, \textrm{aged}}_{\textrm{M},\ell,\alpha} + B^{k, \textrm{aged}}_{\textrm{F},\ell,\alpha}
$$

$$
B^\textrm{aged}_{\ell,\alpha} = B^\textrm{aged}_{\textrm{M},\ell,\alpha} + B^\textrm{aged}_{\textrm{F},\ell,\alpha} 
$$






(apportion-unaged)=
## Apportioning within unaged fish sampless

For unaged fish, the procedure is more complicated, because the samples are not weighed individually, but are weighed only for each sex separately. Therefore, the length-weight relationship derived from aged fish samples (which are weighed individually) is used to partition the total weight into different length bins. Within each length bin, the age distribution of aged fish samples of that length is further used to partition the weight into different length-age bins.

### Unaged weight proportions of different length bins

The number proportion of unaged fish samples of length $\ell$ with respect to all unaged fish samples is

$$
{r_n}^\textrm{unaged/unaged}_\ell = \frac{ n^\textrm{unaged}_\ell }{ n^\textrm{unaged} }
$$

The number proportion of unaged fish samples of sex $s$ and length $\ell$ with respect to all unaged fish samples is

$$
{r_n}^\textrm{unaged/unaged}_{s,\ell} = \frac{ n^\textrm{unaged}_{s,\ell} }{ n^\textrm{unaged} }
$$

Here the number proportions are expressed as empirical probability mass function derived from the unaged fish samples.

Using the length-weight relationship derived from all aged fish samples (both male and female), the weight proportion of unaged fish samples of length $\ell$ with respect to all unaged fish samples is

$$
{r_w}^\textrm{unaged/unaged}_\ell = \frac{ {r_n}^\textrm{unaged/unaged}_\ell \times W(\ell) }{ \sum_\ell {r_n}^\textrm{unaged/unaged}_\ell \times W(\ell) }
$$


### Unaged weight proportions of different sexes
Separately, the weight proportions of male and female within the unaged fish samples can similarly be calculated by using the number proportion of a given sex and the sex-specific length-weight relationship.

The inferred total weights of male and female unaged fish are

$$
w^\textrm{unaged}_\textrm{M} = \sum_\ell {r_n}^\textrm{unaged/unaged}_{M,\ell} \, n^\textrm{unaged}_\textrm{M} \, W_\textrm{M}(\ell)
$$

and

$$
w^\textrm{unaged}_\textrm{F} = \sum_\ell {r_n}^\textrm{unaged/unaged}_{F,\ell} \, n^\textrm{unaged}_\textrm{F} \, W_\textrm{F}(\ell)
$$

The weight proportions of male and female unaged fish with respect to all unaged fish samples are then

$$
{r_w}^\textrm{unaged/unaged}_\textrm{M} = \frac{ w^\textrm{unaged}_\textrm{M} }{ w^\textrm{unaged}_\textrm{M} + w^\textrm{unaged}_\textrm{F} }
$$

and 

$$
{r_w}^\textrm{unaged/unaged}_\textrm{F} = \frac{ w^\textrm{unaged}_\textrm{M} }{ w^\textrm{unaged}_\textrm{M} + w^\textrm{unaged}_\textrm{F} }
$$

The weight proportions of unaged male and female fish with respect to all fish samples (aged _and_ unaged) are then

$$
{r_w}^\textrm{unaged/all}_\textrm{M} = {r_w}^\textrm{unaged/unaged}_\textrm{M} \, {r_w}^\textrm{unaged/all}
$$

and 

$$
{r_w}^\textrm{unaged/all}_\textrm{F} = {r_w}^\textrm{unaged/unaged}_\textrm{F} \, {r_w}^\textrm{unaged/all}
$$


### Unaged biomass apportioned with sex and length

The biomass of unaged fish of sex $s$, length $\ell$ in transect interval $k$ can then be expressed as

$$
B^{k, \textrm{unaged}}_{s, \ell} = \rho^k_B \, A^k \, {r_w}^\textrm{unaged/all}_s \, {r_w}^\textrm{unaged/unaged}_\ell
$$

and the biomass of unaged fish of length $\ell$ in transect interval $k$ can be expressed as

$$
B^{k, \textrm{unaged}}_\ell = \rho^k_B \, A^k \, {r_w}^\textrm{unaged/all} \, {r_w}^\textrm{unaged/unaged}_\ell
$$

Summing across all transect intervals, the total apportioned biomass of unaged fish of sex $s$, length $\ell$ is then

$$
B^\textrm{unaged}_{s, \ell} = \sum_k B^{k, \textrm{unaged}}_{s, \ell} 
$$

and the total apportioned biomass of unaged fish of length $\ell$ is

$$
B^\textrm{unaged}_\ell = \sum_k B^{k, \textrm{unaged}}_\ell
$$


### Unaged biomass apportioned with sex, length, and age

To further partition the biomass of different length bins into different age bins, the weight proportion of fish samples of age $\alpha$ with respect to all aged fish samples of length $\ell$ is used. The calculation is carried out for male and female fish separately before the estimates are combined.

The total apportioned biomass of unaged male and female fish of length $\ell$ and age $\alpha$ is

$$
B^\textrm{unaged}_{\textrm{M},\ell,\alpha} = B^\textrm{unaged}_{\textrm{M},\ell} \times \frac{ B^\textrm{aged}_{\textrm{M},\ell,\alpha} }{ \sum_\alpha B^\textrm{aged}_{\textrm{M},\ell,\alpha} }
$$

$$
B^\textrm{unaged}_{\textrm{F},\ell,\alpha} = B^\textrm{unaged}_{\textrm{F},\ell} \times \frac{ B^\textrm{aged}_{\textrm{F},\ell,\alpha} }{ \sum_\alpha B^\textrm{aged}_{\textrm{F},\ell,\alpha} }
$$

The total apportioned biomass of unaged fish of length $\ell$ and age $\alpha$ is then

$$
B^\textrm{unaged}_{\ell,\alpha} = B^\textrm{unaged}_{\textrm{M},\ell,\alpha} + B^\textrm{unaged}_{\textrm{F},\ell,\alpha}
$$



## Total biomass apportioned with length and age

Using the quantities from above, the total apportioned biomass of all (aged _and_ unaged) fish of length $\ell$ and age $\alpha$ is then

$$
B^\textrm{all}_{\ell,\alpha} = B^\textrm{aged}_{\ell,\alpha} + B^\textrm{unaged}_{\ell,\alpha}
$$
