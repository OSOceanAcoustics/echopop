(apportion-abundance)=
# Apportioning back-calculated abundance

```{attention} 
Back-calculating and apportioning abundance from kriged biomass estimates as described here has not been implemented in `Echopop`.
```


## Back-calculating abundance from kriged biomass estimates

The biomass estimates for male and female fish ($s=M$ and $s=F$, respectively) along transect interval $k$ across all lengths ($\ell$) and all ages ($\alpha$) are:

$$ 
B_{\textrm{M}}^{k} =
\sum_{\ell} B_{\textrm{M}, \ell}^{k, \textrm{unaged}} +
\sum_{\ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, \textrm{aged}}
\label{eq:biomass_M} \tag{1}
$$

$$
B_{\textrm{F}}^{k} =
\sum_{\ell} B_{\textrm{F}, \ell}^{k, \textrm{unaged}} +
\sum_{\ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, \textrm{aged}}
\label{eq:biomass_F} \tag{2}
$$

The biomass estimates for all fish including both sexed and unsexed fish in the transect interval $k$ is then:

$$
B^k = B_\textrm{M}^k + B_\textrm{F}^k.
$$


The estimated abundance $\hat{N}^k$ can be back-calculated from the kriged biomass estimates using an averaged length-weight relationship $\overline{W}(\ell)$ via:

$$
\hat{N}^{k} = \frac{B^{k}}{\overline{W}(\ell)},
\label{eq:abundance} \tag{3}
$$

where $\overline{W}(\ell)$ is the length-weight regression relationship derived from the catch data. 

Similarly, $\hat{\textit{NASC}^{k}}$ can be back-calculated from the estimated abundance using the averaged differential backscattering cross-section of the $i^{\text{th}}$ stratum, $\bar{\sigma}_{bs}^i$, via:

$$
\hat{\textit{NASC}^k} = \hat{N}^k \times \bar{\sigma}_\textrm{bs}^i,
$$

when the transect interval $k$ falls in stratum $i$. See [](stratification) for more information.


```{note} 
In Chu's Echopro implementation, both $\hat{N}_{s}^{k}$ and $\hat{N}^{k}$ are calculated using a single $\overline{W}(\ell)$ fit from **all** (male, female, and unsexed) fish samples, instead of sex-specific fits.
```



## Apportioning back-calculated abundance

Below, the back-calculated $\hat{N}^k$ $\eqref{eq:abundance}$ is apportioned similarly to the [<b>weight proportions</b>](apportion_biomass.md#unaged-biomass-apportioned-with-sex-length-and-age) across sex, length, and age. 


### Number of fish samples

#### Unaged fish

The numbers of unaged male and female fish of length $\ell$ are:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M},\ell}^{\textrm{unaged}} &= \sum_{j \in J_{\textrm{M},\ell}^{\textrm{unaged}}}n_j \nonumber \\
    n_{\textrm{F},\ell}^{\textrm{unaged}} &= \sum_{j \in J_{\textrm{F},\ell}^{\textrm{unaged}}}n_j \nonumber
\end{aligned}
\label{eq:total_unaged_sex_length} \tag{4}
\end{equation}
$$

Therefore, the total numbers of male and female unaged fish of length $\ell$ are:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M}}^{\textrm{unaged}} &= \sum_{\ell}n_{\textrm{M},\ell}^{\textrm{unaged}} \nonumber \\
    n_{\textrm{F}}^{\textrm{unaged}} &= \sum_{\ell}n_{\textrm{F},\ell}^{\textrm{unaged}} \nonumber
\end{aligned}
\label{eq:total_unaged_sex} \tag{5}
\end{equation}
$$


#### Aged fish

The numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M},\ell,\alpha}^{\textrm{aged}} &= \sum_{j \in J_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}n_j \nonumber \\
    n_{\textrm{F},\ell,\alpha}^{\textrm{aged}} &= \sum_{j \in J_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}n_j \nonumber
\end{aligned}
\label{eq:total_aged_sex_length_age} \tag{6}
\end{equation}
$$

Therefore, the total numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

$$
\begin{equation}
\begin{aligned}
    n_{\textrm{M}}^{\textrm{aged}} &= \sum_{\ell,\alpha}n_{\textrm{M},\ell,\alpha}^{\textrm{aged}} \nonumber \\
    n_{\textrm{F}}^{\textrm{aged}} &= \sum_{\ell,\alpha}n_{\textrm{F},\ell,\alpha}^{\textrm{aged}} \nonumber \\    
\end{aligned}
\label{eq:total_aged_sex} \tag{7}
\end{equation}
$$



### Number proportions

The sex-specific numbers for unaged $\eqref{eq:total_unaged_sex}$ and aged $\eqref{eq:total_aged_sex}$ fish are then summed to calculate the total number of unaged fish ($n^{\textrm{unaged}}$), aged ($n^{\textrm{aged}}$), and all ($n$) fish:

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

#### Unaged fish

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

#### Aged fish

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




### Apportioning abundances

#### Unaged fish

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

#### Aged fish

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


#### Combining unaged and aged estimates

Lastly, the estimated abundance of all fish (including unaged and aged fish) of length $\ell$ can be obtained by:

$$
\hat{N}_{\ell}^{k,i} = \hat{N}_{\ell}^{k, \textrm{unaged}} + \sum_{\alpha} \hat{N}_{\ell,\alpha}^{k, \textrm{aged}}
\label{eq:abundance_length} \tag{17}
$$