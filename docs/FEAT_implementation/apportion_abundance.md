(apportion-abundance)=
# Apportioning back-calculated abundance

## Back-calculating abundance from kriged biomass estimates

The biomass estimates for male and female fish ($s=M$ and $s=F$, respectively) along transect interval $k$ across all lengths ($\ell$) and all ages ($\alpha$) are:

(eq-biomass-male)=
$$ 
    B_{\textrm{i,M}}^{k} =
        \sum_{\ell} B_{i,\textrm{M}, \ell}^{k, \textrm{unaged}} +
            \sum_{\ell, \alpha} B_{i,\textrm{M}, \ell, \alpha}^{k, \textrm{aged}}
$$

(eq-biomass-female)=
$$
    B_{i,\textrm{F}}^{k} =
        \sum_{\ell} B_{i,\textrm{F}, \ell}^{k, \textrm{unaged}} +
            \sum_{\ell, \alpha} B_{i,\textrm{F}, \ell, \alpha}^{k, \textrm{aged}}
$$

The biomass estimates for all fish including both sexed and unsexed fish in the transect interval $k$ is then:

(eq-biomass-total)=
$$
    B^k_i = B_{i,\textrm{M}}^k + B_{i,\textrm{F}}^k.
$$


The estimated abundance $\hat{N}^k$ can be back-calculated from the kriged biomass estimates using an averaged length-weight relationship $\overline{\mathcal{W}}(\ell)$ via:

(eq-abundance-total)=
$$
    \hat{N}^{k}_i = \frac{B^{k}_i}{\overline{\mathcal{W}}(\ell)_\text{all}},
$$

where $\overline{W}(\ell)_\text{all}$ is the length-weight regression relationship derived from the catch data. 

Similarly, $\hat{\textit{NASC}_i^{\,\,k}}$ can be back-calculated from the estimated abundance using the averaged differential backscattering cross-section of the $i^{\text{th}}$ stratum, $\bar{\sigma}_{bs}^i$, via:

(eq-nasc)=
$$
    \hat{\textit{NASC}^{\,\,k}}_i = \hat{N}^k_i \times \bar{\sigma}_\textrm{bs}^i,
$$

when the transect interval $k$ falls in stratum $i$.


```{note} 
In Chu's Echopro implementation, both $\hat{N}_{i,s}^{k}$ and $\hat{N}^{k}_i$ are calculated using a single $\overline{\mathcal{W}}(\ell)_\text{all}$ fit from **all** (male, female, and unsexed) fish samples, instead of sex-specific fits.
```
## Apportioning back-calculated abundance

Below, the back-calculated $\hat{N}^k_i$ {ref}`Eq. (3) <eq-abundance-total>` is apportioned similarly to the <b>{ref}`weight proportions <unaged-biomass-apportionment>`</b>.

### Number of fish samples

#### Unaged fish

The numbers of unaged male and female fish of length $\ell$ are:

(eq-unaged-total-sex-length)=
$$
\begin{aligned}
n_{i,\mathrm{M},\ell}^{\mathrm{unaged}} &= \sum_{j\in J_{i,\mathrm{M},\ell}^{\mathrm{unaged}}} n_j \\
n_{i,\mathrm{F},\ell}^{\mathrm{unaged}} &= \sum_{j\in J_{i,\mathrm{F},\ell}^{\mathrm{unaged}}} n_j
\end{aligned}
$$

Therefore, the total numbers of male and female unaged fish of length $\ell$ are:

(eq-unaged-sex-length)=
$$
\begin{aligned}
n_{i,\mathrm{M}}^{\mathrm{unaged}} &= \sum_{\ell} n_{i,\mathrm{M},\ell}^{\mathrm{unaged}} \\
n_{i,\mathrm{F}}^{\mathrm{unaged}} &= \sum_{\ell} n_{i,\mathrm{F},\ell}^{\mathrm{unaged}}
\end{aligned}
$$

#### Aged fish

The numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

(eq-aged-sex-length-age-total)=
$$
\begin{aligned}
n_{i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}} &= \sum_{j\in J_{i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}}} n_j \\
n_{i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}} &= \sum_{j\in J_{i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}}} n_j
\end{aligned}
$$

Therefore, the total numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

(eq-aged-sex-length-total)=
$$
\begin{aligned}
n_{i,\mathrm{M}}^{\mathrm{aged}} &= \sum_{\ell,\alpha} n_{i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}} \\
n_{i,\mathrm{F}}^{\mathrm{aged}} &= \sum_{\ell,\alpha} n_{i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}}
\end{aligned}
$$

### Number proportions

The sex-specific numbers for unaged {ref}`Eq. (6) <eq-unaged-sex-length>` and aged {ref}`Eq. (8) <eq-aged-sex-length-total>` fish are then summed to calculate the total number of unaged fish ($n^{\textrm{unaged}}$), aged ($n^{\textrm{aged}}$), and all ($n$) fish:

(eq-total-counts)=
$$
\begin{aligned}
n_i^{\mathrm{unaged}} &= n_{i,\mathrm{M}}^{\mathrm{unaged}} + n_{i,\mathrm{F}}^{\mathrm{unaged}} \\
n_i^{\mathrm{aged}}   &= n_{i,\mathrm{M}}^{\mathrm{aged}}   + n_{i,\mathrm{F}}^{\mathrm{aged}} \\
n_i                    &= n_i^{\mathrm{unaged}} + n_i^{\mathrm{aged}}
\end{aligned}
$$

#### Unaged fish

The number proportions of male and female unaged fish of length $\ell$ {ref}`Eq. (5) <eq-unaged-total-sex-length>` relative to the sex-specific totals of unaged fish {ref}`Eq. (8) <eq-aged-sex-length-total>` are:

(eq-unaged-number-proportions-sex-length)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged/unaged}} &= \frac{n_{i,\mathrm{M},\ell}^{\mathrm{unaged}}}{n_{i,\mathrm{M}}^{\mathrm{unaged}}} \\
r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged/unaged}} &= \frac{n_{i,\mathrm{F},\ell}^{\mathrm{unaged}}}{n_{i,\mathrm{F}}^{\mathrm{unaged}}}
\end{aligned}
$$

The number proportions of male and female unaged fish of length $\ell$ relative to the total number of fish {ref}`Eq. (9) <eq-total-counts>` are:

(eq-unaged-number-proportions-sex)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged/all}} &= \frac{n_{i,\mathrm{M},\ell}^{\mathrm{unaged}}}{n_i} \\
r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged/all}} &= \frac{n_{i,\mathrm{F},\ell}^{\mathrm{unaged}}}{n_i}
\end{aligned}
$$

The within-unaged-group proportion integrated over sex $s$ in stratum $i$ is:

$$
r_{n,i,\ell}^{\mathrm{unaged/unaged}} = \sum_{s} r_{N,i,s,\ell}^{\mathrm{unaged/unaged}}
$$

The number proportions of male and female unaged fish of length $\ell$ with respect to the total number of fish (unaged and aged combined) are:

(eq-unaged-number-proportions)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged/all}} &= r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged/unaged}}\times r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged/all}} \\
r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged/all}} &= r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged/unaged}}\times r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged/all}}
\end{aligned}
$$

#### Aged fish

Similar to the above, the number of male and female aged fish of length $\ell$ and age $\alpha$ {ref}`Eq. (7) <eq-aged-sex-length-age-total>` relative to the sex-specific totals of aged fish {ref}`Eq. (8) <eq-aged-sex-length-total>` are:

(eq-number-proportions-aged-sex-length-age)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged/aged}} &= \frac{n_{i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}}}{n_{i,\mathrm{M}}^{\mathrm{aged}}} \\
r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged/aged}} &= \frac{n_{i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}}}{n_{i,\mathrm{F}}^{\mathrm{aged}}}
\end{aligned}
$$

The number proportions of male and female aged fish of length $\ell$ and age $\alpha$ relative to the total number of fish {ref}`Eq. (9) <eq-total-counts>` are:

(eq-number-proportions-aged-sex)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged/all}} &= \frac{n_{i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}}}{n_i} \\
r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged/all}} &= \frac{n_{i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}}}{n_i}
\end{aligned}
$$

The number proportions of male and female unaged fish of length $\ell$ and age $\alpha$ with respect to the total number of fish (unaged and aged combined) are:

(eq-number-proportions-aged)=
$$
\begin{aligned}
r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged/all}} &= r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged/aged}}\times r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged/all}} \\
r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged/all}} &= r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged/aged}}\times r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged/all}}
\end{aligned}
$$

### Apportioning abundances

#### Unaged fish

For each transect interval $k$, the total estimated abundance of male, female, and all unaged fish of length $\ell$ are apportioned according to the number proportions in {ref}`Eq. (10) <eq-unaged-number-proportions-sex-length>`:

(eq-abundance-unaged)=
$$
\begin{aligned}
\hat{N}_{i,\mathrm{M},\ell}^{k,\mathrm{unaged}} &= \hat{N}_i^{k}\times r_{N,i,\mathrm{M},\ell}^{\mathrm{unaged}} \\
\hat{N}_{i,\mathrm{F},\ell}^{k,\mathrm{unaged}} &= \hat{N}_i^{k}\times r_{N,i,\mathrm{F},\ell}^{\mathrm{unaged}} \\
\hat{N}_{i,\ell}^{k,\mathrm{unaged}} &= \hat{N}_{i,\mathrm{M},\ell}^{k,\mathrm{unaged}} + \hat{N}_{i,\mathrm{F},\ell}^{k,\mathrm{unaged}}
\end{aligned}
$$

#### Aged fish

Similarly, for each transect interval $k$, the total estimated abundance of male, female, and all aged fish of length $\ell$ and age $\alpha$ are apportioned according to the number proportions in {ref}`Eq. (13) <eq-number-proportions-aged-sex-length-age>`: 

(eq-abundance-aged)=
$$
\begin{aligned}
\hat{N}_{i,\mathrm{M},\ell,\alpha}^{k,\mathrm{aged}} &= \hat{N}_i^{k}\times r_{N,i,\mathrm{M},\ell,\alpha}^{\mathrm{aged}} \\
\hat{N}_{i,\mathrm{F},\ell,\alpha}^{k,\mathrm{aged}} &= \hat{N}_i^{k}\times r_{N,i,\mathrm{F},\ell,\alpha}^{\mathrm{aged}} \\
\hat{N}_{i,\ell,\alpha}^{k,\mathrm{aged}} &= \hat{N}_{i,\mathrm{M},\ell,\alpha}^{k,\mathrm{aged}} + \hat{N}_{i,\mathrm{F},\ell,\alpha}^{k,\mathrm{aged}}
\end{aligned}
$$


#### Combining unaged and aged estimates

Lastly, the estimated abundance of all fish (including unaged and aged fish) of length $\ell$ can be obtained by:

$$
\hat{N}_{\ell}^{k,i} = \hat{N}_{\ell}^{k,\mathrm{unaged},i} + \sum_{\alpha} \hat{N}_{\ell,\alpha}^{k,\mathrm{aged},i}
$$
