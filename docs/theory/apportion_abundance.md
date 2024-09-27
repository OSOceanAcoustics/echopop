(apportion-abundance)=
# Back-calculating abundances from kriged biomass estimates

```{attention} 
`Echopop` currently does support kriged abundance back-calculation from biomass estimates. Refer to the <b>[kriged biomass apportionment](apportion_biomass.md)</b> for more information on how `Echopop` incorporates kriging into population estimates.
```

Biomass estimates for each $s$ ($\textrm{M}$ and $\textrm{F}$) along transect interval $k$ within stratum $i$ are summed across $\ell$ and $\alpha$ via:

$$ 
B_{\textrm{M}}^{k, i} = \sum_{\textrm{M}, \ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, i, \textrm{aged}} +  \sum_{\textrm{M}, \ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, i, \textrm{unaged}}
\label{eq:biomass_M} \tag{1}
$$

$$
B_{\textrm{F}}^{k, i} = \sum_{\textrm{F}, \ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, i, \textrm{aged}} + \sum_{\textrm{F}, \ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, i, \textrm{unaged}}
\label{eq:biomass_F} \tag{2}
$$

Similarly, biomass estimates for all fish ($B^{k, i}$), which is inclusive of both sexed and unsexed fish, are also summed (see: <b>[kriged biomass summation for more details](apportion_biomass.md#total-biomass-apportioned-with-length-and-age)</b>).


These kriged biomass estimates are then converted to sexed ($\hat{N}_{s}^{k, i}$) and total ($\hat{N}^{k, i}$) abundance by using averaged length-weight regression output ($\overline{W}^{~i}(\ell)$) for stratum $i$. Consequently, $\overline{W}^i(\ell)$ can be defined either by using the average length-weight relationship produced from each $i$ or parameterizing $\overline{W}^i(\ell)$ with the mean length ($\bar{\ell}^i$) within each $i$. It is important to note, however, that both $\hat{N}_{s}^{k, i}$ and $\hat{N}^{k, i}$ are calculated using $\overline{W}^{~i}(\ell)$ fit from <b>all</b> individuals (i.e. male, female, and unsexed).

$$
\hat{N}^{k, i} = \frac{B^{k, i}}{\overline{W}^{~i}(\ell)}
\label{eq:abundance} \tag{3}
$$

```{note} 
With $\hat{N}_{\textrm{All}}^{k, i}$ calculated, $\hat{\textit{NASC}^{k, i}}$ can then be back-calculated by using the averaged $i^{\text{th}}$ differential backscattering cross-section ($\bar{\sigma}_{\textrm{bs}}^{i}$):
$
\hat{\textit{NASC}^{k, i}} = \hat{N}_{\textrm{All}}^{k, i} \bar{\sigma}_{\textrm{bs}}^{i}
$
```

## Apportioning the back-calculated abundance estimates

### Summing fish counts

The back-calculated $\hat{N}^{k, i}$ $\eqref{eq:abundance}$ is subsequently apportioned similarly to the [<b>weight proportions</b>](apportion_biomass.md#unaged-biomass-apportioned-with-sex-length-and-age) across sex, length, and age. 

#### Unaged fish

This process is first done across $\ell$ for unaged fish, and both $\ell$ *and* $\alpha$ for aged fish. First, the total number counts for unaged fish ($n_{s,\ell}^{i,~\textrm{unaged}}$):

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M},\ell}^{i,~\textrm{unaged}} &= \sum_{j \in J_{\textrm{M},\ell}^{i,~\textrm{unaged}}}n_j \nonumber \\
    n_{\textrm{F},\ell}^{i,~\textrm{unaged}} &= \sum_{j \in J_{\textrm{F},\ell}^{i,~\textrm{unaged}}}n_j \nonumber
\end{aligned}
\label{eq:total_unaged_sex_length} \tag{4}
\end{equation}
$$

These are then summed across all fish of $s$:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M}}^{i,~\textrm{unaged}} &= \sum_{\textrm{M},\ell}n_{\textrm{M},\ell}^{i,~\textrm{unaged}} \nonumber \\
    n_{\textrm{F}}^{i,~\textrm{unaged}} &= \sum_{\textrm{F}, \ell}n_{\textrm{F},\ell}^{i,~\textrm{unaged}} \nonumber
\end{aligned}
\label{eq:total_unaged_sex} \tag{5}
\end{equation}
$$

#### Aged fish

The total counts for aged fish across both $\ell$ and $\alpha$ are similarly calculated via:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M},\ell,\alpha}^{i,~\textrm{aged}} &= \sum_{j \in J_{\textrm{M},\ell,\alpha}^{i,~\textrm{aged}}}n_j \nonumber \\
    n_{\textrm{F},\ell,\alpha}^{i,~\textrm{aged}} &= \sum_{j \in J_{\textrm{F},\ell,\alpha}^{i,~\textrm{aged}}}n_j \nonumber
\end{aligned}
\label{eq:total_aged_sex_length_age} \tag{6}
\end{equation}
$$

These are then summed across all fish of $s$:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M}}^{i,~\textrm{aged}} &= \sum_{\textrm{M},\ell,\alpha}n_{\textrm{M},\ell,\alpha}^{i,~\textrm{aged}} \nonumber \\
    n_{\textrm{F}}^{i,~\textrm{aged}} &= \sum_{\textrm{F},\ell,\alpha}n_{\textrm{F},\ell,\alpha}^{i,~\textrm{aged}} \nonumber \\    
\end{aligned}
\label{eq:total_aged_sex} \tag{7}
\end{equation}
$$

### Number proportions

The sex-specific abundances for unaged $\eqref{eq:total_unaged_sex}$ and aged $\eqref{eq:total_aged_sex}$ fish are then summed together to calculate the total unaged ($n^{i,~\textrm{unaged}}$), aged ($n^{i,~\textrm{aged}}$), and all ($n^{i}$) fish:

$$
\begin{equation}
\begin{aligned}
    n^{i,~\textrm{unaged}} &= n_{\textrm{M}}^{i,~\textrm{unaged}} + n_{\textrm{F}}^{i,~\textrm{unaged}} \nonumber \\
    n^{i,~\textrm{aged}} &= n_{\textrm{M}}^{i,~\textrm{aged}} + n_{\textrm{F}}^{i,~\textrm{aged}} \nonumber \\
    n^{i} &= n^{i,~\textrm{unaged}} + n^{i,~\textrm{aged}} \nonumber
\end{aligned}
\label{eq:total_counts} \tag{8}
\end{equation}
$$

#### Unaged fish

The number counts of unaged fish within $i$ across $\ell$ for each $s$ $\eqref{eq:total_unaged_sex_length}$ relative to the sex-specific totals $\eqref{eq:total_unaged_sex}$, $r_{n,s,\ell}^{i,~\textrm{unaged/unaged}}$, are:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M},\ell}^{i,~\textrm{unaged/unaged}} &= \frac{n_{\textrm{M},\ell}^{i,~\textrm{unaged}}}{n_{\textrm{M}}^{i,~\textrm{unaged}}} \nonumber \\
    r_{n,~\textrm{F},\ell}^{i,~\textrm{unaged/unaged}} &= \frac{n_{\textrm{F},\ell}^{i,~\textrm{unaged}}}{n_{\textrm{F}}^{i,~\textrm{unaged}}} \nonumber    
\end{aligned}
\label{eq:number_proportions_unaged_sex_length} \tag{9}
\end{equation}
$$

In a similar manner, the unaged fish number counts relative to the sum of unaged and aged number counts $\eqref{eq:total_counts}$, $r_{n,s,\ell}^{i,~\textrm{unaged/all}}$, are then calculated via:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M},\ell}^{i,~\textrm{unaged/all}} = \frac{n_{\textrm{M},\ell}^{i,~\textrm{unaged}}}{n^{i}} \nonumber \\
    r_{n,~\textrm{F},\ell}^{i,~\textrm{unaged/all}} = \frac{n_{\textrm{F},\ell}^{i,~\textrm{unaged}}}{n^{i}} \nonumber
\end{aligned}
\label{eq:number_proportions_unaged_sex} \tag{10}
\end{equation}
$$

The number proportions referencing unaged $\eqref{eq:number_proportions_unaged_sex_length}$ and all $\eqref{eq:number_proportions_unaged_sex}$ fish are then combined to calculate the overall sex-specific number proportions:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M}}^{i,~\textrm{unaged}} &= r_{n,~\textrm{M},\ell}^{i,~\textrm{unaged/unaged}} r_{n,~\textrm{M},\ell}^{i,~\textrm{unaged/all}} \nonumber \\
    r_{n,~\textrm{F}}^{i,~\textrm{unaged}} &= r_{n,~\textrm{F},\ell}^{i,~\textrm{unaged/unaged}} r_{n,~\textrm{F},\ell}^{i,~\textrm{unaged/all}} \nonumber
\end{aligned}
\label{eq:number_proportions_unaged} \tag{11}
\end{equation}
$$

#### Aged fish

Similar to unaged fish, the number counts of aged fish within $i$ across $\ell$ and $\alpha$ for each $s$ $\eqref{eq:total_aged_sex_length_age}$ relative to the sex-specific totals $\eqref{eq:total_aged_sex}$, $r_{n,s,\ell,\alpha}^{i,~\textrm{aged/aged}}$, are:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M},\ell,\alpha}^{i,~\textrm{aged/aged}} &= \frac{n_{\textrm{M},\ell,\alpha}^{i,~\textrm{aged}}}{n_{\textrm{M}}^{i,~\textrm{aged}}} \nonumber \\
    r_{n,~\textrm{F},\ell,\alpha}^{i,~\textrm{aged/aged}} &= \frac{n_{\textrm{F},\ell,\alpha}^{i,~\textrm{aged}}}{n_{\textrm{F}}^{i,~\textrm{aged}}} \nonumber
\end{aligned}
\label{eq:number_proportions_aged_sex_length_age} \tag{12}
\end{equation}
$$

In a similar manner, the unaged fish number counts relative to the sum of unaged and aged number counts $\eqref{eq:total_counts}$, $r_{n,s,\ell,\alpha}^{i,~\textrm{aged/all}}$, are then calculated via:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M},\ell,\alpha}^{i,~\textrm{aged/all}} = \frac{n_{\textrm{M},\ell,\alpha}^{i,~\textrm{aged}}}{n^{i}} \nonumber \\
    r_{n,~\textrm{F},\ell,\alpha}^{i,~\textrm{aged/all}} = \frac{n_{\textrm{F},\ell,\alpha}^{i,~\textrm{aged}}}{n^{i}} \nonumber
\end{aligned}
\label{eq:number_proportions_aged_sex} \tag{13}
\end{equation}
$$

The number proportions referencing unaged $\eqref{eq:number_proportions_aged_sex_length_age}$ and all $\eqref{eq:number_proportions_aged_sex}$ fish are then combined to calculate the overall sex-specific number proportions:

$$
\begin{equation}
\begin{aligned}
    r_{n,~\textrm{M}}^{i,~\textrm{aged}} &= r_{n,~\textrm{M},\ell,\alpha}^{i,~\textrm{aged/aged}} r_{n,~\textrm{M},\ell,\alpha}^{i,~\textrm{aged/all}} \nonumber \\
    r_{n,~\textrm{F}}^{i,~\textrm{aged}} &= r_{n,~\textrm{F},\ell,\alpha}^{i,~\textrm{aged/aged}} r_{n,~\textrm{F},\ell,\alpha}^{i,~\textrm{aged/all}} \nonumber
\end{aligned}
\label{eq:number_proportions_aged} \tag{14}
\end{equation}
$$

### Apportioning abundances

#### Unaged fish

Total unaged fish abundance estimates for $k$ within $i$ $\eqref{eq:abundance}$ are then apportioned for each $s$ across $\ell$, $\hat{N}_{s,\ell}^{k, i,~\textrm{unaged}}$ using the computed number proportions $\eqref{eq:number_proportions_unaged_sex_length}$. The sexed estimates are then summed to compute the total unaged fish abundance estimates, $\hat{N}_{\ell}^{k, i,~\textrm{unaged}}$:

$$
\begin{equation}
\begin{aligned}
    \hat{N}_{\textrm{M},\ell}^{k, i,~\textrm{unaged}} &= \hat{N}^{k, i} r_{n,~\textrm{M}}^{i,~\textrm{unaged}} \nonumber \\
    \hat{N}_{\textrm{F},\ell}^{k, i,~\textrm{unaged}} &= \hat{N}^{k, i} r_{n,~\textrm{F}}^{i,~\textrm{unaged}} \nonumber \\
    \hat{N}_{\ell}^{k, i,~\textrm{unaged}} &= \hat{N}_{\textrm{M},\ell}^{k, i,~\textrm{unaged}} + \hat{N}_{\textrm{F},\ell}^{k, i,~\textrm{unaged}} \nonumber
\end{aligned}
\label{eq:abundance_unaged} \tag{15}
\end{equation}
$$

#### Aged fish

Total unaged fish abundance estimates for $k$ within $i$ $\eqref{eq:abundance}$ are then apportioned for each $s$ across $\ell$ and $\alpha$, $\hat{N}_{s,\ell,\alpha}^{k, i,~\textrm{aged}}$ using the computed number proportions $\eqref{eq:number_proportions_aged_sex_length_age}$. The sexed estimates are then summed to compute the total unaged fish abundance estimates, $\hat{N}_{\ell,\alpha}^{k, i,~\textrm{aged}}$:

$$
\begin{equation}
\begin{aligned}
    \hat{N}_{\textrm{M},\ell,\alpha}^{k, i,~\textrm{aged}} &= \hat{N}^{k, i} r_{n,~\textrm{M}}^{i,~\textrm{aged}} \nonumber \\
    \hat{N}_{\textrm{F},\ell,\alpha}^{k, i,~\textrm{aged}} &= \hat{N}^{k, i} r_{n,~\textrm{F}}^{i,~\textrm{aged}} \nonumber \\
    \hat{N}_{\ell,\alpha}^{k, i,~\textrm{aged}} &= \hat{N}_{\textrm{M},\ell,\alpha}^{k, i,~\textrm{aged}} + \hat{N}_{\textrm{F},\ell,\alpha}^{k, i,~\textrm{aged}} \nonumber
\end{aligned}
\label{eq:abundance_aged} \tag{16}
\end{equation}
$$

#### Combining unaged and aged estimates

Lastly, unaged $\eqref{eq:abundance_unaged}$ and aged $\eqref{eq:abundance_aged}$ abundance estimates can be consolidated to apportion the total abundances across $\ell$ irrespective of $\alpha$:

$$
\hat{N}_{\ell}^{k,i} = \hat{N}_{\ell}^{k, i,~\textrm{unaged}} + \sum_{\alpha} \hat{N}_{\ell,\alpha}^{k, i,~\textrm{aged}}
\label{eq:abundance_length} \tag{17}
$$