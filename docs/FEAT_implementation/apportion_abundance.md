(apportion-abundance)=
# Apportioning back-calculated abundance

## Back-calculating abundance from kriged biomass estimates

The biomass estimates for male and female fish ($s=M$ and $s=F$, respectively) along transect interval $k$ across all lengths ($\ell$) and all ages ($\alpha$) are:

(eq-biomass-male)=
$$ 
    B_{\textrm{M}}^{k} =
        \sum_{\ell} B_{\textrm{M}, \ell}^{k, \textrm{unaged}} +
            \sum_{\ell, \alpha} B_{\textrm{M}, \ell, \alpha}^{k, \textrm{aged}}
    \tag{1}
$$

(eq-biomass-female)=
$$
    B_{\textrm{F}}^{k} =
        \sum_{\ell} B_{\textrm{F}, \ell}^{k, \textrm{unaged}} +
            \sum_{\ell, \alpha} B_{\textrm{F}, \ell, \alpha}^{k, \textrm{aged}}
    \tag{2}
$$

The biomass estimates for all fish including both sexed and unsexed fish in the transect interval $k$ is then:

(eq-biomass-total)=
$$
    B^k = B_\textrm{M}^k + B_\textrm{F}^k.
    \tag{3}
$$


The estimated abundance $\hat{N}^k$ can be back-calculated from the kriged biomass estimates using an averaged length-weight relationship $\overline{W}(\ell)$ via:

(eq-abundance-total)=
$$
    \hat{N}^{k} = \frac{B^{k}}{\overline{W}(\ell)},
    \tag{3}
$$

where $\overline{W}(\ell)$ is the length-weight regression relationship derived from the catch data. 

Similarly, $\hat{\textit{NASC}^{k}}$ can be back-calculated from the estimated abundance using the averaged differential backscattering cross-section of the $i^{\text{th}}$ stratum, $\bar{\sigma}_{bs}^i$, via:

(eq-nasc)=
$$
    \hat{\textit{NASC}^k} = \hat{N}^k \times \bar{\sigma}_\textrm{bs}^i,
    \tag{4}
$$

when the transect interval $k$ falls in stratum $i$. See [](stratification) for more information.


```{note} 
In Chu's Echopro implementation, both $\hat{N}_{s}^{k}$ and $\hat{N}^{k}$ are calculated using a single $\overline{W}(\ell)$ fit from **all** (male, female, and unsexed) fish samples, instead of sex-specific fits.
```



## Apportioning back-calculated abundance

Below, the back-calculated $\hat{N}^k$ {ref}`Eq. (3) <eq-abundance-total>` is apportioned similarly to the <b>{ref}`weight proportions <unaged-biomass-apportionment>`</b>.


### Number of fish samples

#### Unaged fish

The numbers of unaged male and female fish of length $\ell$ are:

(eq-unaged-total-sex-length)=
$$
    \begin{equation}
        \begin{aligned}
            n_{\textrm{M},\ell}^{\textrm{unaged}} &= \sum_{j \in J_{\textrm{M},\ell}^{\textrm{unaged}}}n_j \nonumber \\
            n_{\textrm{F},\ell}^{\textrm{unaged}} &= \sum_{j \in J_{\textrm{F},\ell}^{\textrm{unaged}}}n_j \nonumber
        \end{aligned}
        \tag{5}
    \end{equation}
$$

Therefore, the total numbers of male and female unaged fish of length $\ell$ are:

(eq-unaged-sex-length)=
$$
    \begin{equation}
    \begin{aligned}
        n_{\textrm{M}}^{\textrm{unaged}} &= \sum_{\ell}n_{\textrm{M},\ell}^{\textrm{unaged}} \nonumber \\
        n_{\textrm{F}}^{\textrm{unaged}} &= \sum_{\ell}n_{\textrm{F},\ell}^{\textrm{unaged}} \nonumber
    \end{aligned}
    \tag{6}
    \end{equation}
$$

#### Aged fish

The numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

(eq-aged-sex-length-age-total)=
$$
    \begin{equation}
        \begin{aligned}
            n_{\textrm{M},\ell,\alpha}^{\textrm{aged}} &= \sum_{j \in J_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}n_j \nonumber \\
            n_{\textrm{F},\ell,\alpha}^{\textrm{aged}} &= \sum_{j \in J_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}n_j \nonumber
        \end{aligned}
        \tag{7}
    \end{equation}
$$

Therefore, the total numbers of male and female aged fish of length $\ell$ and age $\alpha$ are:

(eq-aged-sex-length-total)=
$$
    \begin{equation}
        \begin{aligned}
            n_{\textrm{M}}^{\textrm{aged}} &= \sum_{\ell,\alpha}n_{\textrm{M},\ell,\alpha}^{\textrm{aged}} \nonumber \\
            n_{\textrm{F}}^{\textrm{aged}} &= \sum_{\ell,\alpha}n_{\textrm{F},\ell,\alpha}^{\textrm{aged}} \nonumber \\    
        \end{aligned}
        \tag{8}
    \end{equation}
$$

### Number proportions

The sex-specific numbers for unaged {ref}`Eq. (6) <eq-unaged-sex-length>` and aged {ref}`Eq. (8) <eq-aged-sex-length-total>` fish are then summed to calculate the total number of unaged fish ($n^{\textrm{unaged}}$), aged ($n^{\textrm{aged}}$), and all ($n$) fish:

(eq-total-counts)=
$$
    \begin{equation}
        \begin{aligned}
            n^{\textrm{unaged}} &= n_{\textrm{M}}^{\textrm{unaged}} + n_{\textrm{F}}^{\textrm{unaged}} \nonumber \\
            n^{\textrm{aged}} &= n_{\textrm{M}}^{\textrm{aged}} + n_{\textrm{F}}^{\textrm{aged}} \nonumber \\
            n &= n^{\textrm{unaged}} + n^{\textrm{aged}} \nonumber
        \end{aligned}
        \tag{9}
    \end{equation}
$$

#### Unaged fish

The number proportions of male and female unaged fish of length $\ell$ {ref}`Eq. (5) <eq-unaged-total-sex-length>` relative to the sex-specific totals of unaged fish {ref}`Eq. (8) <eq-aged-sex-length-total>` are:

(eq-unaged-number-proportions-sex-length)=
$$
    \begin{equation}
    \begin{aligned}
        {r_N}_{\textrm{M},\ell}^{\textrm{unaged/unaged}} &= \frac{n_{\textrm{M},\ell}^{\textrm{unaged}}}{n_{\textrm{M}}^{\textrm{unaged}}} \nonumber \\
        {r_N}_{\textrm{F},\ell}^{\textrm{unaged/unaged}} &= \frac{n_{\textrm{F},\ell}^{\textrm{unaged}}}{n_{\textrm{F}}^{\textrm{unaged}}} \nonumber    
    \end{aligned}
    \tag{10}
    \end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ relative to the total number of fish {ref}`Eq. (9) <eq-total-counts>` are:

(eq-unaged-number-proportions-sex)=
$$
    \begin{equation}
        \begin{aligned}
            {r_N}_{\textrm{M},\ell}^{\textrm{unaged/all}} = \frac{n_{\textrm{M},\ell}^{\textrm{unaged}}}{n} \nonumber \\
            {r_N}_{\textrm{F},\ell}^{\textrm{unaged/all}} = \frac{n_{\textrm{F},\ell}^{\textrm{unaged}}}{n} \nonumber
        \end{aligned}
        \tag{11}
    \end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ with respect to the total number of fish (unaged and aged combined) are:

(eq-unaged-number-proportions)=
$$
    \begin{equation}
        \begin{aligned}
            {r_N}_{\textrm{M}, \ell}^{\textrm{unaged}/\textrm{all}} &= {r_N}_{\textrm{M},\ell}^{\textrm{unaged/unaged}} \times {r_N}_{\textrm{M},\ell}^{\textrm{unaged/all}} \nonumber \\
            {r_N}_{\textrm{F}, \ell}^{\textrm{unaged}/\textrm{all}} &= {r_N}_{\textrm{F},\ell}^{\textrm{unaged/unaged}} \times {r_N}_{\textrm{F},\ell}^{\textrm{unaged/all}} \nonumber
        \end{aligned}
        \tag{12}
    \end{equation}
$$

#### Aged fish

Similar to the above, the number of male and female aged fish of length $\ell$ and age $\alpha$ {ref}`Eq. (7) <eq-aged-sex-length-age-total>` relative to the sex-specific totals of aged fish {ref}`Eq. (8) <eq-aged-sex-length-total>` are:

(eq-number-proportions-aged-sex-length-age)=
$$
    \begin{equation}
        \begin{aligned}
            {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/aged}} &= \frac{n_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}{n_{\textrm{M}}^{\textrm{aged}}} \nonumber \\
            {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/aged}} &= \frac{n_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}{n_{\textrm{F}}^{\textrm{aged}}} \nonumber
        \end{aligned}
        \tag{13}
    \end{equation}
$$

The number proportions of male and female aged fish of length $\ell$ and age $\alpha$ relative to the total number of fish {ref}`Eq. (9) <eq-total-counts>` are:

(eq-number-proportions-aged-sex)=
$$
    \begin{equation}
        \begin{aligned}
            {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/all}} = \frac{n_{\textrm{M},\ell,\alpha}^{\textrm{aged}}}{n} \nonumber \\
            {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/all}} = \frac{n_{\textrm{F},\ell,\alpha}^{\textrm{aged}}}{n} \nonumber
        \end{aligned}
        \tag{14}
    \end{equation}
$$

The number proportions of male and female unaged fish of length $\ell$ and age $\alpha$ with respect to the total number of fish (unaged and aged combined) are:

(eq-number-proportions-aged)=
$$
    \begin{equation}
        \begin{aligned}
            {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged}/\textrm{all}} &= {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/aged}} \times {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged/all}} \nonumber \\
            {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged}.\textrm{all}} &= r_{n,~\textrm{F},\ell,\alpha}^{\textrm{aged/aged}} \times {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged/all}} \nonumber
        \end{aligned}
    \tag{15}
    \end{equation}
$$

### Apportioning abundances

#### Unaged fish

For each transect interval $k$, the total estimated abundance of male, female, and all unaged fish of length $\ell$ are apportioned according to the number proportions in {ref}`Eq. (10) <eq-unaged-number-proportions-sex-length>`:

(eq-abundance-unaged)=
$$
    \begin{equation}
        \begin{aligned}
            \hat{N}_{\textrm{M},\ell}^{k, \textrm{unaged}} &= \hat{N}^{k} \times {r_N}_{\textrm{M},\ell}^{\textrm{unaged}} \nonumber \\
            \hat{N}_{\textrm{F},\ell}^{k, \textrm{unaged}} &= \hat{N}^{k} \times {r_N}_{\textrm{F},\ell}^{\textrm{unaged}} \nonumber \\
            \hat{N}_{\ell}^{k, \textrm{unaged}} &= \hat{N}_{\textrm{M},\ell}^{k, \textrm{unaged}} + \hat{N}_{\textrm{F},\ell}^{k, \textrm{unaged}} \nonumber \\
        \end{aligned}
    \tag{16}
    \end{equation}
$$

#### Aged fish

Similarly, for each transect interval $k$, the total estimated abundance of male, female, and all aged fish of length $\ell$ and age $\alpha$ are apportioned according to the number proportions in {ref}`Eq. (13) <eq-number-proportions-aged-sex-length-age>`: 

(eq-abundance-aged)=
$$
    \begin{equation}
        \begin{aligned}
            \hat{N}_{\textrm{M},\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}^{k} \times {r_N}_{\textrm{M},\ell,\alpha}^{\textrm{aged}} \nonumber \\
            \hat{N}_{\textrm{F},\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}^{k} \times {r_N}_{\textrm{F},\ell,\alpha}^{\textrm{aged}} \nonumber \\
            \hat{N}_{\ell,\alpha}^{k, \textrm{aged}} &= \hat{N}_{\textrm{M},\ell,\alpha}^{k, \textrm{aged}} + \hat{N}_{\textrm{F},\ell,\alpha}^{k, \textrm{aged}} \nonumber \\
        \end{aligned}
    \tag{17}
    \end{equation}
$$


#### Combining unaged and aged estimates

Lastly, the estimated abundance of all fish (including unaged and aged fish) of length $\ell$ can be obtained by:

$$
    \hat{N}_{\ell}^{k,i} = \hat{N}_{\ell}^{k, \textrm{unaged}} + \sum_{\alpha} \hat{N}_{\ell,\alpha}^{k, \textrm{aged}}
    \tag{18}
$$