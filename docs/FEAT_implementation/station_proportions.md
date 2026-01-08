(station-proportions)=

# Number and weight proportions for fish samples

The challenges associated with apportioning the kriged population estimates are to properly combine the various pieces of information to distribute the estimates into different sex, length, and age groups. This is because fish samples obtained from a haul (trawl) are processed at two different stations that report different biometric data:

(apportion-trawl-stations)=
## Station 1 (*unaged*)
Fish specimens are first sorted according to sex. Their respective standard lengths are then grouped by sex and binned into length distributions ({ref}`Eq. 2.2 <eq-22>`). Haul weights correspond to the bulk total of all unaged fish for each trawl. This results in a bifurcated dataset where the length measurements are stored separately from the catch data. 

### Length
| Haul      | Sex          | Length     | Count               |
|-----------|--------------|------------|---------------------|
| 1         | M            | 10 cm      | 6                   |
| 1         | F            | 10 cm      | 17                  |
| 2         | M            | 11 cm      | 0                   |
| 2         | F            | 11 cm      | 22                  |
| ...       | ...          | ...        | ...                 |

### Weight
| Haul  |  Weight    |
|-------|------------|
| 1     | 12 kg      |
| 2     | 132 kg     |
| 3     | 332 kg     |
| 4     | 0 kg       |
| ...   | ...        |

## Station 2 (*aged*)
Fish specimens are first sorted according to sex. Their respective standard lengths and weights are measured at sea, and otolith ages are determined post-survey onshore. So where station 1 represents bulk measurements of fish, station 2 is specimen-specific. 

| Haul      | Sex          | Age        | Length     | Weight     |
|-----------|--------------|------------|------------|------------|
| 1         | M            | 3          | 10 cm      | 0.15 kg    |
| 1         | M            | 4          | 10 cm      | 0.16 kg    |
| 1         | F            | 2          | 11 cm      | 0.18 kg    |
| 2         | M            | 5          | 11 cm      | 0.20 kg    |
| ...       | ...          | ...        | ...        | ...        |

## Combining the aged and unaged datasets

While the binned lengths and ages from aged data can be used to [directly apportion population estimates](apportion-aged), the unaged measurements need to be [further processed](apportion-unaged) such that they share the same dimensions as the aged data. In practice, [hauls are stratified](stratification-intro) based on their length distributions (as determined by a KS-based cluster analysis) or latitudes (INPFC).

## Overall weight proportion of aged and unaged fish samples

### Aged fish samples
For aged fish samples, the data consists of individual fish specimens with sex, length, and age information. Here, $\text{M}$ and $\text{F}$ denote male and female sexes, respectively. The $\text{aged}$ and $\text{unaged}$ superscripts refer to the station-specific datasets mentioned above. The total weight of aged fish for stratum $i$ is calculated by summing over hauls in that stratum, then over male and female specimens across all length and age bins for each haul, where $g$ is a haul belonging to set $G_i$ (hauls in stratum $i$), and $j$ is an individual specimen belonging to set $J$ (all animals):

$$
\begin{split}
w^\textrm{aged}_i
&= \sum_{g \in G^\textrm{aged}_i} w^\textrm{aged}_{g,i} \\
&= \sum_{g \in G^\textrm{aged}_i} \left( w^\textrm{aged}_{g,\textrm{M},i} + w^\textrm{aged}_{g,\textrm{F},i} \right) \\
&= \sum_{g \in G^\textrm{aged}_i} \left( \sum_{\ell,\alpha} w^\textrm{aged}_{g,\textrm{M},\ell,\alpha,i} + \sum_{\ell,\alpha} w^\textrm{aged}_{g,\textrm{F},\ell,\alpha,i} \right) \\
&= \sum_{g \in G^\textrm{aged}_i} \sum_{\ell,\alpha} \sum_{j \in J^\textrm{aged}_{g,\textrm{M},\ell,\alpha,i}} w_j + \sum_{g \in G^\textrm{aged}_i} \sum_{\ell,\alpha} \sum_{j \in J^\textrm{aged}_{g,\textrm{F},\ell,\alpha,i}} w_j~.
\end{split}
$$

### Unaged fish samples
For unaged fish samples, specimen-specific weights were not measured and are instead represented by their haul catch weights. Therefore, the total summed weight of unaged fish in stratum $i$ is:

$$
w^\textrm{unaged}_i = \sum_{g \in G^\textrm{unaged}_i} \sum_{j \in J^{\textrm{unaged}}_{i,g}} w_j.
$$ 

### Overall weight proportion
The overall weight proportion of aged and unaged fish for stratum $i$ is calculated for both the aged (${r_w}^{\text{aged}/\text{all}}_i$) and unaged (${r_w}^{\text{unaged}/\text{all}}_i$) fish via:

$$
\begin{align*}
    {r_w}^\textrm{aged/all}_i &= \frac{ w^\textrm{aged}_i }{ w^\textrm{aged}_i + w^\textrm{unaged}_i } \\
    {r_w}^\textrm{unaged/all}_i &= \frac{ w^\textrm{unaged}_i }{ w^\textrm{aged}_i + w^\textrm{unaged}_i }.
\end{align*}
$$