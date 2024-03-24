# Apportioning kriged outputs

The fish samples from a trawl (haul) are processed at two stations. At station 1 the sex, length, and weight of the fish samples are measured in a more crude sense. At station 2 the sex, length, weight, and age of each fish are measured, and additional tissue samples are collected.

The challenges associated with apportioning the kriged biomass is to properly combine the various pieces of information to distribute the biomass into different sex, length, and age groups. Importantly, the sex, length, and age distributions of fish vary across strata, which needs to be accounted for in the apportioning procedure.

The apportioning involves different steps for the aged and unaged fish samples, which are described separately in [](apportion-aged) and [](apportion-unaged).


## Symbols and notation

The following symbols and notation will be used in the sections below:

### Transect-related parameters
- $i$ denotes a specific stratum
- $h$ denotes a specific haul
- $k$ denotes a specific transect interval, which belongs to a particular stratum depending on the stratification scheme

### Biometric parameters
- $s$ is the sex of the fish sample and can take values $\textrm{M}$ (male) or $\textrm{F}$ (female)
- $\ell$ is the length of the fish sample
- $\alpha$ is the age of the fish sample
- $J$ denotes the set of fish samples:
    - The superscript describes whether the samples are from all fish samples ($\textrm{all}$), the aged fish samples ($\textrm{aged}$), or unaged fish samples ($\textrm{unaged}$) fish samples. If $h$ is added, it indicates the haul the samples are from.
    - The subscript describes the biometric parameters
    - For example:
        - $J^{\textrm{all}, h}$ is the set of all fish samples from haul $h$
        - $J^\textrm{aged}_{s,\ell,\alpha}$ is the set of aged fish samples of sex $s$, length $\ell$, and age $\alpha$
- $W^i_s(\ell)$ describes the length-weight relationship for a population of fish of sex $s$ in stratum $i$. For example:
    - $W^i_\textrm{male}(\ell)$ is the length-weight relationship for male fish in stratum $i$
    - $W^i(\ell)$ is the length-weight relationship for all fish in stratum $i$

### Weight and weight proportion
- $w_j$ is the weight of fish sample $j$
- $r$ denotes the weight proportion
    - The superscript describes the sample population and the reference population, both of which can be the aged fish samples ($\textrm{aged}$), the unaged fish samples ($\textrm{unaged}$), or all fish samples ($\textrm{all}$)
    - The subscript describes the biometric parameters of the samples
    - For example:
        - $r^\textrm{aged/all}$ is the weight proportion of aged fish with reference to all fish samples
        - $r^\textrm{aged/aged}_{s,\ell}$ is the weight propotion of aged fish of sex $s$ and length $\ell$ with reference to all aged fish samples

### Acoustically derived quantities
- $\rho^k_B$ is the biomass density in transect interval $k$
- $A^k$ is the area covered by transect interval $k$
- $B^k$ is the biomass in transect interval $k$

### Others
- The superscript $\textrm{total}$ denotes a summation over all fish samples or all survey intervals. For example:
    - $w^\textrm{total}$ is the total weight of all fish samples
    - $B^\textrm{total}_{s,\ell}$ is the total apportioned biomass of fish of sex $s$ and length $\ell$





## Overall weight proportion of aged and unaged fish samples

### Aged fish samples
The summed weight of aged fish samples is $\sum_{j \in J^\textrm{aged}} w_j$. The data is given in the form of individual fish samples with sex, length, and age information (in `specimen_df`). Therefore, the summed weight is computed by summing over both male and female aged fish samples of all length and age bins:

$$
\begin{split}
w^\textrm{aged}_\textrm{s}
&= \sum_{\ell,\alpha} w^\textrm{aged}_{\textrm{s},\ell,\alpha} \\
&= \sum_{\ell,\alpha} \sum_{j \in J^\textrm{aged}_{\textrm{s},\ell,\alpha}} w_j
\end{split}
$$
where $s$ can be $\textrm{M}$ or $\textrm{F}$.

The total summed weight of aged fish is then

$$
w^\textrm{aged} = w^\textrm{aged}_\textrm{M} + w^\textrm{aged}_\textrm{F}
$$

### Unaged fish samples
The summed weight of unaged fish samples is $\sum_{j \in J^\textrm{unaged}} w_j$. The data is given in the form of haul catch weight (in `haul_df`). Therefore, the total summed weight of unaged fish is

$$
w^\textrm{unaged} = \sum_{j \in J^{\textrm{unaged},h}} w_j
$$

### Weight proportion
The overall weight propotion of aged and unaged fish are then

$$
r^\textrm{aged/total} = \frac{ w^\textrm{aged} }{ w^\textrm{aged} + w^\textrm{unaged} } 
$$

and

$$
r^\textrm{unaged/total} = \frac{ w^\textrm{unaged} }{ w^\textrm{aged} + w^\textrm{unaged} } 
$$

respectively.





(apportion-aged)=
## Apportioning within aged fish samples

The weight proportions of aged fish samples of sex $s$, length $\ell$, and age $\alpha$ with respect to all aged _and_ unaged fish samples are used to apportion the aged component of the total biomass.

The weight propotions are:

$$
r^\textrm{aged/all}_{s,\ell,\alpha} = \frac{ \sum_{j \in J^\textrm{aged}_{s,\ell,\alpha}} w_j }{ \sum_{j \in J^\textrm{all}} w_j }
$$

The apportioned biomass for transect interval $k$ is

$$
B^{k, \textrm{aged}}_{s,\ell,\alpha} = \rho^k_B A^k r^\textrm{aged/all}_{s,\ell,\alpha}
$$

The apportioned biomass over all transect intervals is then

$$
B^\textrm{total,aged}_{s,\ell,\alpha} = \sum_k B^{k, \textrm{aged}}_{s,\ell,\alpha}
$$

It follows naturally that the combined biomass quantities for male and female fish are:

$$
B^{k, \textrm{aged}}_{\ell,\alpha} = B^{k, \textrm{aged}}_{\textrm{M},\ell,\alpha} + B^{k, \textrm{aged}}_{\textrm{F},\ell,\alpha}
$$

$$
B^\textrm{total,aged}_{\ell,\alpha} = B^\textrm{total,aged}_{\textrm{M},\ell,\alpha} + B^\textrm{total,aged}_{\textrm{F},\ell,\alpha} 
$$




(apportion-unaged)=
## Apportioning within unaged fish sampless