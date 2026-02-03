(adaptive-search)=

# Neighborhood selection and adaptive search

The number of samples $n$ appearing in the predictor of {ref}`Eq. (3.19) <eq-319>` is not fixed globally, but is determined locally for each prediction location $\mathbf{u}$ through an adaptive neighborhood search. Although ordinary kriging places no theoretical restriction on the number of observations used in the estimator, the computational cost grows with $O(n^3)$, and distant samples contribute little information while potentially degrading numerical conditioning. For these reasons, local kriging is typically performed using only samples in the vicinity of each prediction location. The method described here defines this local set in an adaptive but spatially uniform manner so that predictions remain stable in both dense and sparse regions.

Let $x_i \in \mathbb{R}^d$ denote a prediction location and let $\{x_j\}_{j=1}^{N}$ denote the observed sample locations. Distances are defined using the Euclidean metric:

$$
    d_{ij} = \|x_i - x_j\|.
    \tag{3.33}
$$

Let $d_{i(1)} \le d_{i(2)} \le \dots \le d_{i(N)}$ denote these distances sorted in ascending order, with corresponding reordered indices $j_{i(k)}$. Only the first $k_{\max}$ nearest samples are considered as candidates, which bounds computational cost independently of the total dataset size. The candidate set is therefore:

$$
    \mathcal{C}_i = \{ j_{i(1)}, \dots, j_{i(k_{\max})} \}.
    \tag{3.34}
$$

A search radius $r$ defines the region of primary spatial influence. The number of candidates lying inside this radius is:

$$
    n_i = \sum_{k=1}^{k_{\max}} \mathbf{1}(d_{i(k)} < r),
    \tag{3.35}
$$

where $\mathbf{1}(\cdot)$ is the indicator function. This count measures the amount of local support available around $x_i$.

The strategy then adapts the neighborhood size using two user defined parameters: a minimum required number of neighbors $k_{\min}$ and the maximum candidate limit $k_{\max}$. The selected neighborhood $\mathcal{N}_i$ is constructed so that it remains as local as possible while guaranteeing sufficient information for a stable kriging solve. When at least $k_{\min}$ samples lie inside the search radius, the neighborhood is composed exclusively of those local observations:

$$
    \mathcal{N}_i = \{ j_{i(k)} : d_{i(k)} < r \}.
    \tag{3.36}
$$

In this regime, the estimator behaves like a classical fixed radius search. The neighborhood reflects only the local spatial structure, and no distant samples are introduced. When fewer than $k_{\min}$ samples fall within the radius but at least one local sample exists, the neighborhood is augmented with the closest exterior samples until the minimum requirement is satisfied. Formally, the neighborhood becomes:

$$
    \mathcal{N}_i = \{ j_{i(1)}, \dots, j_{i(k_{\min})} \}.
    \tag{3.37}
$$

This preserves all available local information while borrowing the smallest possible amount of additional support. The prediction therefore degrades smoothly from interpolation toward mild extrapolation as sampling density decreases, rather than failing abruptly. In extremely sparse regions where no samples lie inside the radius, the method reduces to a pure $k$ nearest neighbor search,

$$
    \mathcal{N}_i = \{ j_{i(1)}, \dots, j_{i(k_{\min})} \}.
    \tag{3.38}
$$

Although these samples may be distant, enforcing the minimum count ensures that the kriging matrix remains well posed and that weights can still be computed reliably. Once $\mathcal{N}_i$ has been determined, the kriging system is assembled in the usual way using only those indices. The weights $\{\lambda_j\}$ and subsequent prediction $\hat{Z}(x_i)$ are identtical to the standard formulation used in ordinary kriging. 

Conceptually, this approach combines the advantages of both fixed radius and fixed $k$ nearest strategies. In dense areas it behaves like a radius based method and uses strictly local information. In sparse areas it behaves like a nearest neighbor method and guarantees numerical stability. Because the neighborhood size adapts automatically to sampling density, predictions remain spatially consistent without requiring manual tuning across different regions of the domain. For this reason, the procedure should be viewed as an implementation policy for neighborhood selection rather than as part of the kriging theory itself. The statistical model and optimality properties of the estimator remain unchanged. The search strategy simply provides a practical and robust way to choose the subset of observations that define each local solve.
