# Apportioning biological variables

<!-- ```{image} ../images/length_stratification.jpg
:alt: length - stratification
:width: 200px
```

```{image} ../images/length_sex_stratification.jpg
:alt: length - sex - stratification
:width: 200px
```

```{image} ../images/length_age_stratification.jpg
:alt: length - age - stratification
:width: 200px
```

```{image} ../images/length_age_sex_stratification.jpg
:alt: length - age- sex - stratification
:width: 200px
``` -->


Both weight ($w_{\alpha, \ell, s}^{i}$) and individual counts ($n_{\alpha, \ell, s}^{i}$) are apportioned in different ways to account for differences among sexes, ages, length, and strata. This also includes slightly different calculations depending on station ($\vartheta$) where animals were either unaged (Station 1, $\vartheta = 1$) or aged (Staion 2, $\vartheta = 2$). In the case of Station 1 where animals are unaged:

$$
n_{\ell, s, \vartheta=1}^{i} = \sum\limits_{n^{i}} x_{\ell, s, \vartheta=1}^{i}
$$

where $x_{\ell, s}^{i}$ represents an individual of sex $s$ belonging to stratum $i$ that has a length within the length bin $\ell$. Once binned, the modeled weight based on the fitted length-weight regression ($\hat{W}_{\ell}$) calculated from all individuals is used to convert $n_{\ell, s, \vartheta=1}^{i}$ into ($\hat{w}_{\ell, s, \vartheta=1}^{i}$) via:

$$
\hat{w}_{\ell, s, \vartheta=1}^{i} = n_{\ell, s, \vartheta=1}^{i} \hat{W}_{\ell}
$$

This is similarly done for $n_{\alpha , \ell, s, \vartheta=2}^{i}$ and $\hat{w}_{\ell, s, \vartheta=2}^{i}$ via:

$$
n_{\alpha, \ell, s, \vartheta=2}^{i} = \sum\limits_{n^{i}} x_{\alpha, \ell, s, \vartheta=2}^{i}
$$

$$
\hat{w}_{\alpha, \ell, s, \vartheta=2}^{i} = n_{\alpha, \ell, s, \vartheta=2}^{i} \hat{W}_{\ell}
$$

In both cases, these can be convered to either number/count or weight proportions for each stratum. For numbered proportions ($p_{n}$):

$$
p_{n: \alpha, \ell, s, \vartheta}^{i} = 
\begin{cases}
    \frac{ n_{\ell, s, \vartheta}^{i} }
         { \sum\limits_{n_{\vartheta}^{i}} n_{\ell, s, \vartheta}^{i} } &
         \text{if } \vartheta = 1 \\
    \frac{ n_{\alpha , \ell, s, \vartheta}^{i} }
         { \sum\limits_{n_{\vartheta}^{i}} n_{\alpha , \ell, s, \vartheta}^{i} } &
         \text{if } \vartheta = 2
\end{cases}
$$

Similarly for the weight proportions ($p_{\hat{w}}$):

$$
p_{\hat{w}: \alpha, \ell, s, \vartheta}^{i} = 
\begin{cases}
    \frac{ \hat{w}_{\ell, s, \vartheta}^{i} }
         { \sum\limits_{w_{\vartheta}^{i}} \hat{w}_{\ell, s, \vartheta}^{i} } &
         \text{if } \vartheta = 1 \\
    \frac{ \hat{w}_{\alpha , \ell, s, \vartheta}^{i} }
         { \sum\limits_{\hat{w}_{\vartheta}^{i}} \hat{w}_{\alpha , \ell, s, \vartheta}^{i} } &
         \text{if } \vartheta = 2
\end{cases}
$$

In some cases, both summed counts, lengths, and proportions from both stations will need to be combined, which requires
adjusting how the proportions are calculated to account for both aged and unaged animals. This first requires summing the total
weights of animals captured in each haul ($\nu_{\vartheta}^{h,i}$) to provide appropriate apportionment:

$$
\begin{align}
    \nu_{\vartheta}^{i} &= \sum\limits_{n_{\vartheta}^{h,i}} \nu_{\vartheta}^{h,i} \\
    \nu^{i} &= \sum\limits_{n_{\vartheta}^{i}}  \nu_{\vartheta}^{i}
\end{align}
$$

For the case where $\vartheta = 2$ where animals are aged, $p_{\hat{w}: \alpha, \ell, s, \vartheta}^{i}$ is
normalized to the summed haul weights (indexed by strata):

$$
\hat{p}_{\hat{w}: \alpha, \ell, s, \vartheta}^{i} =
    \hat{p}_{\hat{w}: \alpha, \ell, s, \vartheta}^{i}
    \left\{
        \frac{ \sum\limits_{n_{\vartheta=2}^{i}} \hat{w}_{\alpha, \ell, s, \vartheta=2}^{i} }
             { \nu^{i} }        
    \right\}
$$

The proportions of aged ($\hat{\vartheta}=\mathrm{T}$) and unaged ($\hat{\vartheta}=\mathrm{F}$) animals are calculated via:

$$
\begin{align}
    p_{\hat{\vartheta}=\mathrm{T}}^{i} &=
        \frac{ \sum\limits_{n_{\vartheta=2}^{i}} \hat{w}_{\alpha, \ell, s, \vartheta=2}^{i} }
            { \nu^{i} } \\
    p_{\hat{\vartheta}=\mathrm{F}}^{i} &=
        1.0 - p_{\hat{\vartheta}=\mathrm{T}}^{i}
\end{align}
$$

For station 1, the length-weight regression for each sex ($\hat{W}_{\ell, s}$) is used to interpolate weights by sex per haul and stratum
($\hat{w}_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{h,i}$):

$$
\hat{w}_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{h,i}: L_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{h,i} \rightarrow \hat{W}_{\ell, s}
$$

This is then summed across strata to get the summed sex-specific stratum weight:

$$
\hat{w}_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{i} =
    \sum\limits_{n_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{i}} \hat{w}_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{h,i}
$$

Once these are calculated, then the weight proportions for $\vartheta=1$ can be calculated. This first requires normalizing the weight
measurements for eahc stratum in station 1:

$$
\hat{\nu}_{\vartheta=1}^{i} =
    \frac{ \nu_{\vartheta=1}^{i} \hat{w}_{\ell, s, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{i} }
         { \hat{w}_{\ell, s=\mathrm{male}, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{i} +
            \hat{w}_{\ell, s=\mathrm{female}, \vartheta=1, \hat{\vartheta}=\mathrm{F}}^{i} }
$$

This is then converted as the normalized proportion:

$$
\hat{p}_{w:s,\vartheta=1}^{i} =
    \frac{\hat{\nu}_{s,\vartheta=1}^{i}}
         { \nu^{i} }
$$

The sexed proportion is then calculated:

$$
\hat{p}_{w:\vartheta=1, s \in \{ \mathrm{male}, \mathrm{female} \}}^{i}=
    \frac{\hat{p}_{w:s,\vartheta=1}^{i}}
         { \hat{p}_{w:s\in \{ \mathrm{male}, \mathrm{female} \},\vartheta=1}^{i} }
$$

This is finally weighted by the unaged proprtions:  

$$
\hat{p}_{w: s, \hat{\vartheta}=\mathrm{F}}^{i} =
    \hat{p}_{w:\vartheta=1, s \in \{ \mathrm{male}, \mathrm{female} \}}^{i} p_{\hat{\vartheta}=\mathrm{F}}^{i}    
$$

which provides the proportion of unaged biomass per sex. 

This then allows the kriged biomass estimates ($\hat{B}_{s}^{i,j,k}$) to be apportioned. The kriged
areal biomass densities ($\hat{\rho}_{A,B}^{i,j,k}$) are first weighted with the above calculated proportions:

$$
\begin{align}
    \hat{\rho}_{A,B: s, \hat{\vartheta}=\mathrm{T}}^{i,j,k} &= 
        \hat{\rho}_{A,B}^{i,j,k} \hat{p}_{\hat{w}: \alpha, \ell, s, \vartheta=2}^{i} \\
    \hat{B}_{s, \hat{\vartheta}=\mathrm{T}}^{i,j,k} &=
        \hat{\rho}_{A,B: s}^{i,j,k} A^{i,j,k}
\end{align}
$$

This quantity represents the apportioned aged biomass. This can be similarly expanded to the unaged biomass calculations:

$$
\begin{align}
    \hat{\rho}_{A,B: s, \hat{\vartheta}=\mathrm{F}}^{i,j,k} &= 
        \hat{\rho}_{A,B}^{i,j,k} \hat{p}_{\hat{w}: \alpha, \ell, s, \vartheta=1}^{i} \\
    \hat{B}_{s, \hat{\vartheta}=\mathrm{F}}^{i,j,k} &=
        \hat{\rho}_{A,B: s}^{i,j,k} A^{i,j,k}
\end{align}
$$

The total biomass is thus estimated via:
$$
\hat{B}_{s}^{i,j,k} =
    \hat{B}_{s, \hat{\vartheta}=\mathrm{T}}^{i,j,k} + \hat{B}_{s, \hat{\vartheta}=\mathrm{F}}^{i,j,k}
$$
