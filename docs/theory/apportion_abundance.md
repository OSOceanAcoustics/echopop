(apportion-abundance)=
# Back-calculating and apportioning abundance estimates


```{attention} 
`Echopop` currently does not support back-calculating abundance from kriged biomass estimates detailed in [](apportion-biomass).
```

```{note}
It is worth noting that all calculations are done for each stratum, $i$. Refer to [](stratification) for more information.
```

Biomass estimates for each $s$ ($\textrm{M}$ and $\textrm{F}$) along transect interval $k$ are summed across $\ell$ and $\alpha$ via:

$$ 
B_{\textrm{M}}^{k} = \sum_{\textrm{M}, \ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, \textrm{aged}} +  \sum_{\textrm{M}, \ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, \textrm{unaged}}
\label{eq:biomass_M} \tag{1}
$$

$$
B_{\textrm{F}}^{k} = \sum_{\textrm{F}, \ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, \textrm{aged}} + \sum_{\textrm{F}, \ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, \textrm{unaged}}
\label{eq:biomass_F} \tag{2}
$$

Similarly, biomass estimates for all fish ($B^{k}$), which is inclusive of both sexed and unsexed fish, are also summed, i.e.,
$$
B^k = B_\textrm{M}^k + B_\textrm{F}^k.
$$


These kriged biomass estimates are then converted to sexed ($\hat{N}_{s}^{k}$) and total ($\hat{N}^{k}$) abundance by using an averaged length-weight relationship ($\overline{W}(\ell)$). $\overline{W}(\ell)$ can be defined either by using the length-weight regression relationship or a parameterized relationship based on mean leangth ($\bar{\ell}$) derived from the catch data. It is important to note, however, that both $\hat{N}_{s}^{k}$ and $\hat{N}^{k}$ are calculated using a $\overline{W}(\ell)$ fit from <b>all</b> individuals (i.e. male, female, and unsexed).

$$
\hat{N}^{k} = \frac{B^{k}}{\overline{W}(\ell)}
\label{eq:abundance} \tag{3}
$$

```{note} 
With $\hat{N}_k$, $\hat{\textit{NASC}^{k}}$ can be back-calculated by using the averaged differential backscattering cross-section for the $i^{\text{th}}$ stratum, $\bar{\sigma}_\textrm{bs}^i$, as
$$
\hat{\textit{NASC}^k} = \hat{N}^k \times \bar{\sigma}_\textrm{bs}^i
$$
```

Below, the back-calculated $\hat{N}^k$ $\eqref{eq:abundance}$ is apportioned similarly to the [<b>weight proportions</b>](apportion_biomass.md#unaged-biomass-apportioned-with-sex-length-and-age) across sex, length, and age. 


## Number of fish samples

### Unaged fish

The number of unaged fish of sex $s$ and length $\ell$ is ($n_{s,\ell}^{\textrm{unaged}}$) of length $\ell$ is:

$$
\begin{equation}
    n_{s,\ell}^{\textrm{unaged}} = \sum_{j \in J_{s,\ell}^{\textrm{unaged}}}n_j,
\label{eq:total_unaged_sex_length} \tag{4}
\end{equation}
$$
where $s=M$ and $s=F$ indicate male and female fish, respectively.

Therefore, the total number of fish of sex $s$ is:

$$
\begin{equation}
    n_s^{\textrm{unaged}} = \sum_\ell n_{s,\ell}^{\textrm{unaged}}
\label{eq:total_unaged_sex} \tag{5}
\end{equation}
$$

### Aged fish

The total number of fish of sex $s$, length $\ell$, and age $\alpha$ is similarly:

$$
\begin{equation}
    n_{s,\ell,\alpha}^{\textrm{aged}} = \sum_{j \in J_{s,\ell,\alpha}^{\textrm{aged}}} n_j
\label{eq:total_aged_sex_length_age} \tag{6}
\end{equation}
$$

The number of aged fish of sex $s$ is then:

$$
\begin{equation}
    n_s^{\textrm{aged}} = \sum_{\ell,\alpha}n_{s,\ell,\alpha}^{\textrm{aged}}
\label{eq:total_aged_sex} \tag{7}
\end{equation}
$$



## Number proportions

The sex-specific numbers for unaged $\eqref{eq:total_unaged_sex}$ and aged $\eqref{eq:total_aged_sex}$ fish are then summed to calculate the total number of unaged fish($n^{\textrm{unaged}}$), aged ($n^{\textrm{aged}}$), and all ($n$) fish:

$$
\begin{equation}
\begin{aligned}
    n^{\textrm{unaged}} &= n_{\textrm{M}}^{\textrm{unaged}} + n_{\textrm{F}}^{\textrm{unaged}} \nonumber \\
    n^{\textrm{aged}} &= n_{\textrm{M}}^{\textrm{aged}} + n_{\textrm{F}}^{\textrm{aged}} \nonumber \\
    n &= n^{\textrm{unaged}} + n^{\textrm{aged}} \nonumber
\end{aligned}
\label{eq:total_counts} \tag{8}
\end{equation}
$$

### Unaged fish

The number proportions of male and female unaged fish of length $\ell$ $\eqref{eq:total_unaged_sex_length}$ relative to the sex-specific totals of unaged fish $\eqref{eq:total_unaged_sex}$ are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M},\ell}^{\textrm{unaged/unaged}} &= \frac{n_{\textrm{M},\ell}^{\textrm{unaged}}}{n_{\textrm{M}}^{\textrm{unaged}}} \nonumber \\
    {r_N}_{\textrm{F},\ell}^{\textrm{unaged/unaged}} &= \frac{n_{\textrm{F},\ell}^{\textrm{unaged}}}{n_{\textrm{F}}^{\textrm{unaged}}} \nonumber    
\end{aligned}
\label{eq:number_proportions_unaged_sex_length} \tag{9}
\end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ relative to the total number of fish $\eqref{eq:total_counts}$ are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M},\ell}^{\textrm{unaged/all}} = \frac{n_{\textrm{M},\ell}^{\textrm{unaged}}}{n} \nonumber \\
    {r_N}_{\textrm{F},\ell}^{\textrm{unaged/all}} = \frac{n_{\textrm{F},\ell}^{\textrm{unaged}}}{n} \nonumber
\end{aligned}
\label{eq:number_proportions_unaged_sex} \tag{10}
\end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ with respect to the total number of fish (unaged and aged combined) are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M}, \ell}^{\textrm{unaged}/\textrm{all}} &= {r_N}_{\textrm{M},\ell}^{\textrm{unaged/unaged}} \times {r_N}_{\textrm{M},\ell}^{\textrm{unaged/all}} \nonumber \\
    {r_N}_{\textrm{F}, \ell}^{\textrm{unaged}/\textrm{all}} &= {r_N}_{\textrm{F},\ell}^{\textrm{unaged/unaged}} \times {r_N}_{\textrm{F},\ell}^{\textrm{unaged/all}} \nonumber
\end{aligned}
\label{eq:number_proportions_unaged} \tag{11}
\end{equation}
$$

### Aged fish

Similar to the above, the number of male and female aged fish of length $\ell$ and age $\alpha$ $\eqref{eq:total_aged_sex_length_age}$ relative to the sex-specific totals of aged fish $\eqref{eq:total_aged_sex}$ are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/aged}} &= \frac{n_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}{n_{\textrm{M}}^{\textrm{aged}}} \nonumber \\
    {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/aged}} &= \frac{n_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}{n_{\textrm{F}}^{\textrm{aged}}} \nonumber
\end{aligned}
\label{eq:number_proportions_aged_sex_length_age} \tag{12}
\end{equation}
$$

The number proportions of male and female aged fish of length $\ell$ and age $\alpha$ relative to the total number of fish $\eqref{eq:total_counts}$ are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/all}} = \frac{n_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}{n} \nonumber \\
    {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/all}} = \frac{n_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}{n} \nonumber
\end{aligned}
\label{eq:number_proportions_aged_sex} \tag{13}
\end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ and age $\alpha$ with respect to the total number of fish (unaged and aged combined) are:

$$
\begin{equation}
\begin{aligned}
    {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged}/\textrm{all}} &= {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/aged}} \times {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/all}} \nonumber \\
    {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged}.\textrm{all}} &= r_{n,~\textrm{F},\ell,\alpha}^{\textrm{aged/aged}} \times {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/all}} \nonumber
\end{aligned}
\label{eq:number_proportions_aged} \tag{14}
\end{equation}
$$




## Apportioning abundances

### Unaged fish

For each transect interval $k$, the total estimated abundance of male, female, and all unaged fish of length $\ell$ are apportioned according to the number proportions in $\eqref{eq:number_proportions_unaged_sex_length}$:

$$
\begin{equation}
\begin{aligned}
    \hat{N}_{\textrm{M},\ell}^{k, \textrm{unaged}} &= \hat{N}^{k} \times {r_N}_{\textrm{M},\ell}^{\textrm{unaged}} \nonumber \\
    \hat{N}_{\textrm{F},\ell}^{k, \textrm{unaged}} &= \hat{N}^{k} \times {r_N}_{\textrm{F},\ell}^{\textrm{unaged}} \nonumber \\
    \hat{N}_{\ell}^{k, \textrm{unaged}} &= \hat{N}_{\textrm{M},\ell}^{k, \textrm{unaged}} + \hat{N}_{\textrm{F},\ell}^{k, \textrm{unaged}} \nonumber \\
\end{aligned}
\label{eq:abundance_unaged} \tag{15}
\end{equation}
$$

### Aged fish

Similarly, for each transect interval $k$, the total estimated abundance of male, female, and all aged fish of length $\ell$ and age $\alpha$ are apportioned according to the number proportions in $\eqref{eq:number_proportions_aged_sex_length_age}$: 

$$
\begin{equation}
\begin{aligned}
    \hat{N}_{\textrm{M},\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}^{k} \times {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged}} \nonumber \\
    \hat{N}_{\textrm{F},\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}^{k} \times {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged}} \nonumber \\
    \hat{N}_{\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}_{\textrm{M},\ell,\alpha}^{k, \textrm{aged}} + \hat{N}_{\textrm{F},\ell,\alpha}^{k, \textrm{aged}} \nonumber \\
\end{aligned}
\label{eq:abundance_aged} \tag{16}
\end{equation}
$$


### Combining unaged and aged estimates

Lastly, the estimated abundance of all fish (including unaged and aged fish) of length $\ell$ can be obtained by:

$$
\hat{N}_{\ell}^{k,i} = \hat{N}_{\ell}^{k, \textrm{unaged}} + \sum_{\alpha} \hat{N}_{\ell,\alpha}^{k, \textrm{aged}}
\label{eq:abundance_length} \tag{17}
$$