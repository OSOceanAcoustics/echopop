# Apportioning kriged outputs

The fish samples from a trawl (haul) are processed at two stations. At station 1 the sex, length, and weight of the fish samples are measured in a more crude sense. At station 2 the sex, length, weight, and age of each fish are measured, and additional tissue samples are collected.

The challenges associated with apportioning the kriged biomass is to properly combine the various pieces of information to distribute the biomass into different sex, length, and age groups. Importantly, the sex, length, and age distributions of fish vary across strata, which needs to be accounted for in the apportioning procedure.

The apportioning involves different steps for the aged and unaged fish samples.


## Symbols and notation

The following symbols and notation will be used in the sections below:

- $s$ is the sex of the fish sample and can take values $\textrm{M}$ (male) or $\textrm{F}$ (female)
- $\ell$ is the length of the fish sample
- $\alpha$ is the age of the fish sample
- $J$ denotes the set of fish samples:
    - The superscript describes whether the samples are from the aged ($\textrm{aged}$) or unaged ($\textrm{unaged}$) fish samples, and from which haul ($h$)
    - The subscript describes the biometric parameters
    - For example, $J^\textrm{aged}_{s,\ell,\alpha}$ is the set of aged fish samples of sex $s$, length $\ell$, and age $\alpha$.
- $w_j$ is the weight of fish sample $j$
- $r$ denotes the weight proportion
    - The superscript describes the sample population and the reference population, both of which can be the aged fish samples ($\textrm{aged}$), the unaged fish samples ($\textrm{unaged}$), or all fish samples ($\textrm{all}$)
    - The subscript describes the biometric parameters of the samples
    - For example, $r^\textrm{aged/all}$ is the weight proportion of aged fish with reference to all fish samples, and $r^\textrm{aged/aged}_{s,\ell}$ is the weight propotion of aged fish of sex $s$ and length $\ell$ with reference to all aged fish samples
- $\rho^k_B$ is biomass density in survey interval $k$
- $A^k$ is the area covered by survey interval $k$
- $B$ is biomass
- Superscript $\textrm{total}$ denotes a summation over all fish samples or all survey intervals. For example, $w^\textrm{total}$ is the total weight of all fish samples, and $B^\textrm{total}_{s,\ell}$ is the total apportioned biomass of fish of sex $s$ and length $\ell$.
- $h$ denotes a specific haul
- $i$ denotes a specific stratum


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
The summed weight of unaged fish samples is $\sum_{j \in J^\textrm{unaged}} w_j$. The data is given in the form of haul catch weight (in `length_df`). Therefore, the total summed weight of unaged fish is

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


## Apportioning within aged fish samples
The weight proportions of aged fish samples of sex $s$, length $\ell$, and age $\alpha$ with respect to all aged fish samples are used to apportion the aged component of the total biomass.




## Apportioning within unaged fish sampless