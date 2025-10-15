# Expressions for transect processing steps

## Data binning

```python
utils.binify(...)
```

When given a value for age ($a$), the bin assignment ($\alpha$) is:

$$
a \in (\alpha_{j-1} - \Delta \alpha,~ \alpha_j + \Delta \alpha]
\quad \text{where} \quad (\alpha_{j-1} - \Delta \alpha) < a \leq (\alpha_j + \Delta \alpha)~.
$$

Length values ($L$) are similarly binned ($\ell$):

$$
L \in (\ell_{i-1} - \Delta \ell, ~\ell_i + \Delta \ell] \quad \text{where} \quad (\ell_{i-1} - \Delta \ell) < L \leq (\ell_i + \Delta \ell)~.
$$

Here, $\Delta \alpha$ and $\Delta \ell$ represent the half-width of the bins, calculated as the average difference between consecutive bin edges.
where:

The distributions of $\alpha$ and $\ell$ are represented as $\mathbf{\vec{\alpha}}$ and $\mathbf{\vec{\ell}}$, respectively:

$$
\mathbf{\vec{\ell}} = \begin{bmatrix}
\ell_1 \\
\ell_2 \\
\ell_3 \\
\vdots
\end{bmatrix}~,
$$

$$
\mathbf{\vec{\alpha}} = \begin{bmatrix}
\alpha_1 \\
\alpha_2 \\
\alpha_3 \\
\vdots
\end{bmatrix}~.
$$

## Length-weight regression fitting

```python
biology.fit_length_weight_regression(...)
```

The log-linear relationship between specimen wet weight ($w$) and $L$ is modeled via:

$$
\log_{10}\big(w(L)\big) = \beta_0 + \beta_1 \cdot \log_{10}(L)~.
$$

Here, $\beta_0$ and $\beta_1$ are the regression coefficients (intercept and slope), and the summation minimizes the squared error between the observed and predicted log-transformed weights:

Here, $\beta_0$ and $\beta_1$ are the regression coefficients, representing the intercept and slope of the log-linear relationship between $w$ and $L$. The equation minimizes the sum of squared differences (errors) between the observed log-transformed weights $\log_{10}(w_i)$ and the predicted values $\beta_0 + \beta_1 \log_{10}(L_i)$ across all specimens ($i = 1, \dots, n$):

$$
(\beta_0, \beta_1) = \underset{(\beta_0, \beta_1)}{\argmin}
\sum_{i=1}^{n} \left( \log_{10}(w_i) - (\beta_0 + \beta_1 \log_{10}(L_i)) \right)^2~.
$$

This process can also be performed separately for sex $s$, resulting in sex-specific regression coefficients $\beta_{0,s}$ and $\beta_{1,s}$.

## Quantize counts

```python
get_proportions.compute_binned_counts(...)
```

The number of animals in a given combination of age $\alpha$, length $\ell$, sex $s$, and stratum $h$ is computed as:

$$
n_g = \sum_{i \in g} 1~,
$$

where $n_g$ represents the total count of specimens belonging to the multi-dimensional bin $g$, where $g$ can include combinations of variables such as age ($\alpha$), length ($\ell$), sex ($s$), and stratum ($h$). For instance, aged fish can be distributed across the multi-dimensional bin $(\alpha, \ell, s, h)$ where:

$$
g = (\alpha, \ell, s, h)~,
$$

whereas unaged fish can be distributed across $(\ell, s, h)$ where:

$$
g = (\ell, s, h)~.
$$

Variables included within $g$ are specified via:

```python
get_proportions.compute_binned_counts(..., groupby_cols)
```

So for aged fish, `groupby_cols=["stratum_num", "length_bin", "age_bin", "sex"]` translates to:

$$
n_{\alpha, \ell, s, h} = \sum_{i \in (\alpha, \ell, s, h)} 1~.
$$

Similarly for unaged fish, `groupby_cols=["stratum_num", "length_bin", "sex"]` translates to:

$$
n_{\ell, s, h} = \sum_{i \in (\ell, s, h)} 1~.
$$

## Calculate the mean weight per length bin

```python
biology.length_binned_weights(...)
```

The fitted average weights per length bin for all, $\hat{W}(\ell)$, and sex-specific, $\hat{W}_s(\ell)$, fish can be computed via two different methods: 1) calculating mean $w_\ell$ and $w_{s, \ell}$ and 2) applying the log-linear regression coefficients. The latter is only used when:

```python
biology.length_binned_weights(..., minimum_threshold_count, impute_bins)
```

When `impute_bins=True`, `minimum_threshold_count` (default: `0`) set a threshold for when fitted weights are imputed for particular $\ell$. For example, `minimum_threshold_count=5` would result in the applying the fitted log-linear regression to $\ell$ where the bin counts ($n_\ell$) are less (or equal to in the case of `minimum_threshold_count = 0)`) this minimum threshold. Otherwise, the mean weight is used. This can be expressed by:

$$
\hat{W}(\ell) =
\begin{cases}
10^{\beta_0} \cdot \ell^{\beta_1}, & \text{if } n_{\ell} < \text{threshold}, \\
\frac{1}{n_{\ell}} \sum\limits_{i \in \ell} w_i, & \text{if } n_{\ell} \geq \text{threshold}.
\end{cases}~,
$$

where $\hat{W}(\ell)$ is the fitted average weight for each $\ell$. Since only the aged data are used for this fitting, $n_\ell$ is therefore:

$$
n_\ell = \sum_{\alpha} \sum_{s} n_{\alpha, \ell, s}~,
$$

when computed for all fish.

This calculation is also done for each $s$ where:

$$
\hat{W}_s(\ell) =
\begin{cases}
10^{\beta_{0,s}} \cdot \ell^{\beta_{1,s}}, & \text{if } n_{\ell, s} < \text{threshold}, \\
\frac{1}{n_{\ell, s}} \sum\limits_{i \in \ell, s} w_i, & \text{if } n_{\ell, s} \geq \text{threshold}.
\end{cases}~,
$$

where:

$$
n_{\ell, s} = \sum_{\alpha} n_{\alpha, \ell, s}~.
$$

## Number proportions

```python
get_proportions.number_proportions(...)
```

Number proportions ($\pi$) are calculated both within and across each dataset. These proportions can be expressed as $\pi_g^{c/C}$ where the proportion of counts in each multi-dimensional bin $g$ in sub-category $c$ (e.g. aged fish) is normalized by the total count across the complete set of categories $C$ (e.g. aged or all fish). This is expressed via:

$$
\pi_g^{c/C} = \frac{n_g^c}{\displaystyle \sum_{c \in C} \sum_{g \in c} n_g^{c}}
\quad \text{for } g \in c~.
$$

The category $c$ is defined via:

```python
get_proportions.number_proportions(..., column_aliases)
```

which maps each $c$ to the resulting output dictionary of `get_proportions.number_proportions`. So `column_aliases = ["aged", "unaged"]` would define $c$ as $\text{aged}$ and $\text{uanged}$ for the calculations of $\pi^{\text{aged}/C}$ and $\pi^{\text{unaged}/C}$, respectively.

### Grouped/partitioned number proportions

Proportion calculations can also be partitioned across subsets of the reference set $C$ by specifying one or more grouping variables (i.e. specific column names):

```python
get_proportions.number_proportions(..., group_columns)
```

The previous equations assume the total counts in the denominator are aggregated over the entirety of $C$, without regard to subgroupings. However, grouping variables instruct the function to compute the proportions within each level of the specified grouping. For instance, `group_columns = ["stratum_num"]` would calculate separate proportions for each stratum $h$.

To formalize this, the subset of $C$ relevant for relevant for a given multi-dimensional bin $g$ is defined as:

$$
C_{\gamma}(g) :=
\begin{cases}
\{ g^* \in C : g^*_{\gamma} = g_{\gamma} \}, & \text{if } \gamma \neq \varnothing \\[8pt]
C, & \text{if } \gamma = \varnothing
\end{cases}~,
$$

where $\gamma$ denotes the grouping variable, representing a subset of the dimensions of $g$. For example, if $\gamma = \{\text{stratum}\}$, then $g_\gamma$ is the stratum number ($h$) value of $g$. So if $g_\gamma = h = 7$, this means $g$ belongs to stratum 7.

The subset $C_{\gamma}(g)$ is then defined as all bins $g^* \in C$ whose grouping variable components $g^*_\gamma$ exactly match $g_\gamma$. In the example above, this means $C_\gamma(g) \implies C_{h=7}(g)$, which thus maps to all multi-dimensional bins within $C$ belonging to $h=7$. Therefore, the subset $C_\gamma(g)$ partitions the reference set $C$ into groups based on those defined by `group_columns`. When `group_columns = []` or `group_columns = None`, meaning $\gamma = \varnothing$, $C_\gamma(g)$ defaults to the entire reference set $C$ with no partitioning.

The grouped proportion is then computing by normalizing counts to the total over subset $C_{\gamma}(g)$:

$$
\pi_g^{c/C} =
\frac{n_g^c}
     {\displaystyle \sum_{c \in C_\gamma} \sum_{g^* \in C_{\gamma}(g)} n_{g^*}^c}
\quad \text{for } g \in c~.
$$

### Filtered number proportions
In some cases, specific groups may be intentionally excluded from the output by applying filters to omit their contributions to the numerator of the proportion calculations. This is done via:

```python
get_proportions.number_proportions(..., exclude_filters)
```

The `exclude_filters` argument effectively removes any bins matching the given criteria **from the numerator only**, while leaving the **denominator untouched**. That is, excluded values are omitted from the output proportions but are still included in the total count used for normalization. For example, if `exclude_filters = {"sex": "unsexed"}`, then all bins where $s = \text{unsexed}$ are excluded from the output and do not appear in any numerator $n_g^c$. However, these same bins are still included in the denominator, since they are part of the full dataset used to define $C$, and are therefore still summed across all $g^* \in C_\gamma(g)$.

To formalize this, let $\mathcal{E} \subset C$ represent the set of bins excluded via `exclude_filters`. Then the grouped proportion is redefined only using non-excluded bins, while the denominator remains unchanged. The condition $g \in C \setminus \mathcal{E}$ restricts the output to non-excluded entries, meaning $n_g^c$ appears in the numerator only if $g \notin \mathcal{E}$. However, $n_{g^*}^c$ in the denominator is summed across all bins $g^* \in C_\gamma(g)$. This ensures that the calculated proportions are normalized over the complete dataset, even if filtered components are excluded from the results. Therefore, $\pi_g^{c/C}$ can be refined as:

$$
\pi_g^{c/C} =
\frac{n_g^c}
     {\displaystyle \sum_{c \in C_\gamma} \sum_{g^* \in C_{\gamma}(g)} n_{g^*}^c}
\quad \text{for } g \in c \setminus \mathcal{E}~.
$$

Because excluded bins are omitted from the numerator but still included in the denominator, the resulting number proportions may **not sum to 1.00** across the output bins. This is expected given that the proportions now represent the share of each non-excluded bin relative to the total count (inclusive of the excluded bins) used for normalization.

### Application to the aged/unaged samples

The aged number proportions incorporate a sex-based exclusion filter defined by:

$$
\mathcal{E} := \{ (\alpha, \ell, s, h) \in \text{aged} : s = \text{unsexed}\}~,
$$

which results in the proportion calculations being limited to only $s \in \{\text{female},
\text{male}\}$ for the within- ($\text{aged}/\text{aged}$) and across-group ($\text{aged}/\text{all}$). Both the aged and unaged proportions are also calculated for each stratum where:

$$
\begin{aligned}
\gamma &= h \\
C_\gamma(g) &\implies C_h(g)
\end{aligned}~.
$$

Two different proportions are calculated for the aged and unaged fish: (1) within-group (`proportion`) and (2) across-group (`proportion_overall`). For aged fish, $c=\text{aged}$, and $C=\{\text{aged}\}$ and $C=\{\text{aged, unaged}\}$ for the within- and across-group proportions, respectively. Thus:

$$
\begin{aligned}
\pi^{\text{aged}/\text{aged}}_{\alpha, \ell, s, h} &=
    \frac
    {
        n^{\text{aged}}_{\alpha, \ell, s, h}
    }
    {
        \sum\limits_{g^*\in \{\text{aged}\}_h(\alpha, \ell, s, h)}
        n^{\text{aged}}_{g^*}
    }
    \quad \text{for } (\alpha, \ell, s, h) \in \text{aged} \setminus \mathcal{E} \\
\pi^{\text{aged}/\text{all}}_{\alpha, \ell, s, h} &=
    \frac
    {
        n^{\text{aged}}_{\alpha, \ell, s, h}
    }
    {
        \sum\limits_{{c \in \{\text{aged}, \text{unaged}\}}_h}
        \sum\limits_{g^* \in \{\text{aged}, \text{unaged}\}_h(\ell, s, h)}
        n^{c}_{g^*}
    }
    \quad \text{for } (\alpha, \ell, s, h) \in \text{aged} \setminus \mathcal{E}
\end{aligned}~.
$$

Unlike the aged proportion calculations, the unaged proportions do not incorporate $\mathcal{E}$; however, they are otherwise calculated similarly. For unaged fish, $c=\text{unaged}$, and $C=\{\text{unaged}\}$ and $C=\{\text{aged, unaged}\}$ for the within- and across-group proportions, respectively. Thus:

$$
\begin{aligned}
\pi^{\text{unaged}/\text{unaged}}_{\ell, s, h} &=
    \frac
    {
        n^{\text{unaged}}_{\ell, s, h}
    }
    {
        \displaystyle \sum_{g^* \in \{\text{unaged}\}_h(\ell, s, h)} n^{\text{unaged}}_{g^*}
    } \\
\pi^{\text{unaged}/\text{all}}_{\ell, s, h} &=
    \frac
    {
        n^{\text{unaged}}_{\ell, s, h}
    }
    {
        \sum\limits_{{c \in \{\text{aged}, \text{unaged}\}}_h}
        \sum\limits_{g^* \in \{\text{aged}, \text{unaged}\}_h(\ell, s, h)}
        n^{c}_{g^*}
    }
\end{aligned}~.
$$

## Distribute weights over defined bins

```python
get_proportions.binned_weights(...)
```

Weights can be summed across multi-dimensional bins $g$ via:

$$
w_{g} = \sum_{i \in g} w_i~.
$$

### Interpolating binned weights

The argument `interpolate` indicates whether or not weights are interpolated across $g$ via:

```python
get_proportions.binned_weights(..., length_weight_dataset, interpolate=True)
```

The length-fitted weights, $W(\ell)$ and $W_s(\ell)$, can be used to estimate weights ($\hat{w}$) when weight measurements are unavailable. Since the $W(\ell)$ and $W_s(\ell)$ were fitted based on the center $L$ of each $\ell$, weights can be interpolated for values of $L$:

$$
\hat{w}(L) = \text{Interp}\big(L, W(\ell)\big)~.
$$

These values can then be summed similar to the previous expressions via:

$$
\hat{w}_{g} = \sum_{i \in g} \hat{w}(L_i)~.
$$

### Filtered binned weights
In some cases, this weight bin summation may only be desired for specific data slices and multi-dimensional bins $g$. This type of inclusive filter can be defined via:

```python
get_proportions.binned_weights(..., include_filter)
```

To formalize this, let $\mathcal{F}$ represent the set of specimens that satisfy the inclusion filter. The summation is then performed only over specimens that meet these criteria. The condition $i \in g \cap \mathcal{F}$ restricts the summation to only those specimens $i$ that are both in the target group $g$ and satisfy the inclusion filter. Therefore, the filtered sum of weights $w_g$ is defined as:

$$
w_{g} = \sum_{i \in g \cap \mathcal{F}} w_i~.
$$

For example, if `include_filter = {"sex": ["female", "male"]}`, then only bins where $s \in \{\text{male, female}\}$ are included in the equations and therefore the outputs. This would therefore be expressed via:

$$
\mathcal{F} = {i \mid s_i \in {\text{male, female}}}~.
$$

### Application to the aged/unaged samples
The binned aged weights are calculated via:

$$
w_{\alpha, \ell, s, h} =
\sum_{i \in (\alpha, \ell, s, h) \cap \mathcal{F}} w_i~.
$$

The binned unaged weights are calculated via:

$$
\hat{w}_{\ell, s, h} =
n^{\text{unaged}}_{\ell, s, h}
\sum_{i \in (\ell, s, h) \cap \mathcal{F}} \hat{w}(L_i) ~.
$$

## Stratum-averaged weights

```python
get_proportions.stratum_averaged_weight(...)
```

Stratum-averaged weights are calculated by combining $\pi_{\alpha, \ell, s, h}^{\text{aged}/\text{aged}}$/$\pi_{\alpha, \ell, s, h}^{\text{aged}/\text{all}}$ and $\pi_{\ell, s, h}^{\text{unaged}/\text{unaged}}$/$\pi_{\ell, s, h}^{\text{unaged}/\text{all}}$ with $W(\ell)$. This first involves calculating the relative proportions of each group:

$$
\pi^{c/\text{all}}_{s,h} =
    \sum_{g \in C_h(s,h)}
    \pi^{c/\text{all}}_g~,
$$

where $s \in \{\text{all}, \text{female}, \text{male}\}$ and:

$$
\sum_{c \in C_h(h)} \pi^{c/\text{all}}_{s,h} = 1.0 ~.
$$

These overall proportions are then used to re-weight and standardize the weight proportions for $c \in \{\text{aged}, \text{unaged}\}$. This involves first adjusting the overall unaged proportions:

$$
\hat{\pi}^{\text{unaged}/\text{all}}_{s,h}=
\frac
{\sum_s \pi^{\text{unaged}/\text{all}}_{s, h}}
{\pi^{\text{unaged}/\text{all}}_{s, h} + \sum_s \pi^{\text{unaged}/\text{all}}_{s, h}}~.
$$

The quantity $\hat{\pi}^{\text{unaged}/\text{all}}_{s,h}$ is then used to adjust the overall aged proportions:

$$
\hat{\pi}^{\text{aged}/\text{all}}_{s,h}=
\frac
{\pi^{\text{aged}/\text{all}}_{s, h}}
{\pi^{\text{aged}/\text{all}}_{s, h} + \hat{\pi}^{\text{unaged}/\text{all}}_{s,h}}~.
$$

The within-group proportions for $c=\text{aged}$ and $c=\text{unaged}$ (distributed over $\ell$ for each $s$ and $h$) are then calculated via:

$$
\pi_{\ell, s, h}^{c/c} =
\sum_{g \in C_h(\ell, s, h)} \pi^{c/c}_g~.
$$

The adjusted overall across- and within-group proportions are then combined to compute the average weights within each stratum. When $s = \text{all}$:

$$
\hat{W}_h = \hat{W}(\ell) \cdot
\left[
        \pi_{\ell, s, h}^{\text{aged}/\text{aged}} \times
        \sum_s \pi^{\text{aged}/\text{all}}_{s, h} +
        \pi_{\ell, s, h}^{\text{unaged}/\text{unaged}} \times
        \sum_s \pi^{\text{unaged}/\text{all}}_{s, h}
\right]~.
$$

Otherwise, a similar calculation is done for $s \in \{\text{female}, \text{male}\}$ using the adjusted proportion calculations via:

$$
\hat{W}_{s, h} =
    \hat{W}_s(\ell) \cdot
    \left[
        \pi_{\ell, s, h}^{\text{aged}/\text{aged}} \times
        \sum_s \hat{\pi}^{\text{aged}/\text{all}}_{s,h} +
        \pi_{\ell, s, h}^{\text{unaged}/\text{unaged}} \times
        \sum_s \hat{\pi}^{\text{aged}/\text{all}}_{s,h}
    \right]~.
$$

## (Aged) Weight proportions

```python
get_proportions.weight_proportions(...)
```

The `pandas.DataFrame` input for `catch_data` includes data collected from each net trawl $t$. The column `haul_weight` corresponds to the total unaged weight ($w_t$) for that particular haul. Since these represent only unaged specimens, the sum of each strata is calculated via:

$$
w^{\text{unaged}}_h = \sum\limits_{ t \in \textbf{t}(h) } w_t~,
$$

where $\textbf{t}$ represents the complete set of hauls.

The aged stratum weights are calculated directly from the input `pandas.DataFrame` (`weight_data`) via:

$$
w^{\text{aged}}_h = \sum_{g \in \{\text{aged}\}_h(\alpha, \ell, s, h)} w_g~,
$$

The total stratum weight is therefore defined as:

$$
w^{\text{all}}_h = w^{\text{unaged}}_h + w^{\text{aged}}_h
$$

The aged weight proportions can be then be calculated using:

$$
\omega^{\text{aged}/\text{all}}_{\alpha, \ell, s, h} =
\frac
{
    w_{\alpha, \ell, s, h}
}
{
    w^{\text{all}}_h
}
$$

## (Unaged) Standardized weight sums

```python
get_proportions.standardize_weights_by_stratum(...)
```

The fitted unaged weights, $\hat{w}(L_i)$, are standardized to better correspond to $w^{\text{unaged}}_h$. This operation recomputes the total unaged weights per stratum for each sex via:

$$
\tilde{w}^{\text{unaged}}_{s, h}=
\left(
\frac
{
    \sum\limits_{g \in \{\text{unaged}\}_h(s, h)} \hat{w}_g
}
{
    \sum\limits_{g \in \{\text{unaged}\}_h(h)} \hat{w}_g
}
\right)
w^{\text{unaged}}_h~.
$$


## (Unaged) Standardized weight proportions

```python
get_proportions.standardize_weight_proportions(...)
```

The standardized unaged weight proportions, $\tilde{w}^{\text{unaged}}_{s, h}$, are used to calculate a new total stratum weight via:

$$
\tilde{w}^\text{all}_h = w^\text{aged}_h + \sum\limits_s \tilde{w}^{\text{unaged}}_{s, h}~.
$$

These stratum totals are then used to calculate the overall sexed weight proportions

$$
\tilde{\omega}^{\text{unaged}/\text{all}}_{s,h}=
\frac
{
    \tilde{w}^{\text{unaged}}_{s, h}
}
{
    \tilde{w}^\text{all}_h
}~.
$$

The within-group unaged weight proportions summed for each sex are then back-calculated from $\tilde{\omega}^{\text{unaged}/\text{all}}_{s,h}$:

$$
\tilde{\omega}^{\text{unaged}/\text{unaged}}_{s,h}=
\frac
{
    \tilde{\omega}^{\text{unaged}/\text{all}}_{s,h}
}
{
    \sum\limits_s \tilde{\omega}^{\text{unaged}/\text{all}}_{s,h}
}~.
$$

The within-group weight proportions across length bins are generated using the associated number proportions, $\pi^{\text{unaged}/\text{unaged}}_{\ell, s, h}$ and fitted weights computed from all fish, $\hat{W}(\ell)$. First, the average fitted weight per length bin is computed via:

$$
\hat{w}_{\ell, h}^{\text{unaged}} = \pi^{\text{unaged}/\text{unaged}}_{\ell, s, h} \cdot \hat{W}(\ell)~,
$$

which is in turn used to compute the within-group proportions:

$$
\omega^{\text{unaged}/\text{unaged}}_{\ell, h}=
\frac
{
    \hat{w}_{\ell, h}^{\text{unaged}}
}
{
    \sum\limits_{\ell} \hat{w}_{\ell, h}^{\text{unaged}}
}~.
$$

Next, the overall proportion of unaged weights per stratum can be simply calculated by:

$$
\omega^{\text{unaged}/\text{all}}_{h} =
1.00 - \sum\limits_{\alpha, \ell, s} \omega^{\text{aged}/\text{all}}_{\alpha, \ell, s, h}~.
$$

This enables the calculation of the overall length-binned weight proportions via:

$$
\omega^{\text{unaged}/\text{all}}_{\ell, h}=
\omega^{\text{unaged}/\text{unaged}}_{\ell, h} \times
\omega^{\text{unaged}/\text{all}}_{h}~.
$$

The last component of this calculation finally involves distributed the sex proportions over $\omega^{\text{unaged}/\text{all}}_{\ell, h}$:

$$
\omega^{\text{unaged}/\text{all}}_{\ell, s, h}=
\omega^{\text{unaged}/\text{all}}_{\ell, h} \times
\tilde{\omega}^{\text{unaged}/\text{unaged}}_{s,h}~.
$$
