(kriging-theory)=
# Interpolating spatial fields

## Ordinary kriging

**Ordinary kriging** is a **Best-Linear-Unbiased Predictor (BLUP)** for spatial interpolation that combines observed values with a model of spatial dependence to produce point estimates and associated uncertainties. The method requires two essential stationarity assumptions.[^cressie] 

1. **First-Order Stationarity**: the mean of the field is unknown but constant over the local prediction neighborhood. 
2. [**Intrinsic Stationarity**](semivariogram_theory.md#the-intrinsic-model): the spatial dependence is only a function of lag vector $\mathbf{h}$, not location $\mathbf{x}$.

This form of kriging is preferred for locally stationary fields or when no reliable trend model is apparent. Since ordinary kriging assumes a locally constant but unknown mean, it is mathematically constructed to filter out this mean. When the regional mean trends slowly, or where any measured trend is very gradual or unreliable, then the small search neighborhood becomes small enough such that the local mean appears effectively constant. Other kriging variants exist (for example, methods that model an explicit spatial trend or include external covariates), and they differ primarily in how the mean is represented and estimated. Mentioning these alternatives is useful for methodological context, but the equations below are restricted to ordinary kriging only.[^cressie][^chiles]

--- 

## General assumptions and principles

The underlying random process [$Z(\mathbf{x})$](./semivariogram_theory.md#eq-randomz) can be characterized by both deterministic and stochastic components:

$$
Z(\mathbf{x}) = m(\mathbf{x}) + \varepsilon(\mathbf{x}),
\tag{1a}
\label{eq:kriging1a}
$$

$$
\mathbb{E}[Z(\mathbf{x})] = m(\mathbf{x}),
\tag{1b}
\label{eq:kriging1b}
$$

$$
\mathbb{E}[\varepsilon(\mathbf{x})] = 0,
\tag{1c}
$$

where $\mathbb{E}[\cdot]$ is the expectation, $m$ is the mean at location $\mathbf{x}$, and $\varepsilon$ represents the combined unobserved small-scale processes and observation errors. This represents the **Universal Kriging** model, which can account for clearly defined global or regional trends. However, estimating both the spatial correlation (i.e. variogram of the residuals, $\varepsilon$) *and* the coefficients of the trend function $m(\mathbf{x})$ simultaneously can become very complex and computationally taxing.

:::{admonition} Theoretical vs. empirical expressions
:class: note
While $Z(\mathbf{x})$ denotes the random field, spatially varying values can also be expressed as $z(\mathbf{x})$ when discussing realized observations at location $\mathbf{x}$, or $\mathbf{z}^*(\mathbf{x})$ for predicted realizations at location $\mathbf{x}$. This distinction makes clear whether an expression is a theoretical statement or an empirical computation.
:::

### Simplification for Ordinary Kriging

For Ordinary Kriging, this complex, position-dependent mean $m(\mathbf{x})$ is simplified by assuming that, within the local neighbor used for estimation, the mean is a constant scalar $\mathbf{m}$, although its specific value is unknown. This assumption of first-order stationarity therefore simplifies $Z(\mathbf{x})$ to:

$$
Z(\mathbf{x}) = \mathbf{m} + \varepsilon(\mathbf{x}).
\tag{2a}
$$

The expectation of the field within the local area is now $\mathbf{m}$, thus:

$$
\mathbb{E}[Z(\mathbf{x})] = \mathbf{m}.
\tag{2b}
$$

### Applying the Intrinsic Model

The assumption of local mean stationarity, combined with the intrinsic model concerning the variance of the residuals, defines the spatial correlation structured used in Ordinary Kriging. This results in two important conditions. First, the expectation of the increments of $Z(\mathbf{x})$ is constant and zero:

$$
\mathbb{E}[Z(\mathbf{x}+\mathbf{h}) - Z(\mathbf{x})] = 0.
\tag{3a}
$$

Second, the variance of the increments depends solely on the lag vector $\mathbf{h}$ (distance and direction) and not on the absolute location $\mathbf{x}$:

$$
\mathbb{V}[Z(\mathbf{x}+\mathbf{h}) - Z(\mathbf{x})] = 2\gamma(\mathbf{h}),
\tag{3b}
$$

where $\gamma(\mathbf{h})$ is the semivariogram.

For processes that satisfy the stronger condition of second-order stationarity (i.e. finite variance), the covariance function $\mathbb{C}(\mathbf{h})$ is defined using the expectation of the product of the field values separated by $\mathbf{h}$:

$$
\mathbb{C}(\mathbf{h}) = \mathbb{E}[Z(\mathbf{x})Z(\mathbf{x}+\mathbf{h})]-\mathbf{m}^2.
\tag{4}
$$

### Anisotropy

In many spatial datasets, the correlation structure is direction-dependent, a property known as anisotropy. Ordinary kriging accounts for anisotropy by transforming the coordinate system so that distances are rescaled according to the principal axes of spatial continuity. Specifically, the lag vector $\mathbf{h}$ between two locations is projected onto these axes, yielding components $\Delta x$ and $\Delta y$. The anisotropic distance is then defined as

$$ 
d_{\text{anisotropic}}^2 = \left(\frac{\Delta x}{a_x}\right)^2 + \left(\frac{\Delta y}{a_y}\right)^2, 
\tag{5}
$$

where $\Delta x$ and $\Delta y$ are the components of $\mathbf{h}$ along the principal axes, and $a_x$ and $a_y$ are the corresponding range parameters. The principal axes are determined by fitting an anisotropic variogram model, which may involve rotating the coordinate system so that the $x$-axis aligns with the direction of maximum spatial continuity and the $y$-axis with the direction of minimum continuity. This transformation ensures that the variogram or covariance model accurately reflects the spatial structure of the data, and the kriging system is constructed using these rescaled distances.

---

## Predictor and unbiasedness

:::{admonition} Note on Notation
:class: note
The variable $\mathbf{x}$ is used in Equations (1-4) to denote a general, theoretical location when defining the properties of the random field. From this point forward, the variable $\mathbf{u}$ is used to denote the target location. The indices $\mathbf{u}_b$, $i$, and $j$ all reference the same set of $\mathbf{n}$ sampled data locations ($\mathbf{u}_1, \mathbf{u}_2, \dots, \mathbf{u}_n$). The indices $i$ (rows) and $j$ (columns) are used for indexing the Kriging matrix rows and columns because they correspond to locations $\mathbf{u}_i$ and $\mathbf{u}_j$, respectively. All indices run from 1 to $n$. 
:::

The general form of the linear predictor, using the $\mathbf{u}$ notation for the prediction location and the observed data $z(\mathbf{u}_b)$, includes the unknown local mean $\mathbf{m}$:

$$
\mathbf{z}^*(\mathbf{u}) = \sum_{b=1}^n \lambda_b(\mathbf{u}) z(\mathbf{u}_b) + 
\left[ 
    1 - \sum_{b=1}^n \lambda_b(\mathbf{u}) 
\right] \mathbf{m}.
\label{eq:linearpredictor}
\tag{6}
$$

In this equation, $n$ represents the number of data points ($z(\mathbf{u}_b)$) selected from the entire dataset that are used to calculate the prediction at a single location $\mathbf{u}$. This number $n$ is a **user-defined** parameter determined by a local search strategy (e.g. maximum distance radius, maximum number of neighbors) implemented during the kriging process. For instance, the search neighborhood can be constrained based on a maximum distance radius informed by a multiple (e.g. 1x, 2x, 3x) of the semivariogram correlation range ($a$). The variable $b$ serves as a generic iterator for summing over these $n$ data points, and $\lambda_b(\mathbf{u})$ are the unknown weights associated with each point.

### Unbiasedness and optimality

A major advantage conferred by using Ordinary Kriging lies in the **unbiasedness constraint**. To ensure the prediction is unbiased regardless of $\mathbf{m}$, the coefficient preceding $\mathbf{m}$ in $\eqref{eq:linearpredictor}$ must be zero. This therefore requires:

$$
\sum_{b=1}^{n} \lambda_b(\mathbf{u}) = 1,
\tag{7a}
$$

which results with the because the bracketed term in $\eqref{eq:linearpredictor}$ being zeroed out:

$$
\left[
    1 - \sum_{b=1}^n \lambda_b(\mathbf{u}) 
\right] 
\mathbf{m} = \left[ 
    1 - 1 
\right] \mathbf{m} = 0.
\tag{7b}
$$

Consequently, the unbiasedness constraint eliminates $\mathbf{m}$ from $\mathbf{z}^*(\mathbf{u})$, therefore simplifying $\eqref{eq:linearpredictor}$ to:

$$
\mathbf{z}^*(\mathbf{u}) = \sum_{b=1}^n \lambda_b(\mathbf{u}) z(\mathbf{u}_b). 
\label{eq:kriging_pred}
\tag{8}
$$

A second condition, **optimality**, requires that the weights minimize the estimation variance, $\sigma_{E}^{2}(\mathbf{u})$, among all linear unbiased estimators:

$$
\min \sigma_E^2(\mathbf{u}) = \min \mathbb{V}[\mathbf{z}^*(\mathbf{u}) - Z(\mathbf{u})].
\tag{9}
\label{eq:sigma_min}
$$

To solve this constrained optimization problem, the spatial relationships between all $n$ data locations must be quantified in a matrix, which defines the initial $n \times n$ covariance structure ($\mathbf{\Gamma}_{n \times n}$):

$$
\begin{equation}
\mathbf{\Gamma}_{n \times n} = 
\begin{bmatrix}
    \gamma(\mathbf{u}_1 - \mathbf{u}_1) & \gamma(\mathbf{u}_1 - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_1 - \mathbf{u}_n) \\
    \gamma(\mathbf{u}_2 - \mathbf{u}_1) & \gamma(\mathbf{u}_2 - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_2 - \mathbf{u}_n) \\
    \vdots & \vdots & \ddots & \vdots \\
    \gamma(\mathbf{u}_n - \mathbf{u}_1) & \gamma(\mathbf{u}_n - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_n - \mathbf{u}_n)
\end{bmatrix},
\tag{10}
\label{eq:kriging_matrix}
\end{equation}
$$

where $\gamma(\mathbf{u}_i - \mathbf{u}_j)$ is the semivariance between data points $i$ and $j$. The minimized $\sigma_E^2$ $\eqref{eq:sigma_min}$ is calculated using the weights and values from $\mathbf{\Gamma}_{n \times n}$. The full expression of $\sigma_E^2(\mathbf{u})$ includes a double summation over all weights and their corresponding semivariances to account for the spatial dependence between every pair of data points, $\mathbf{u}_j$ and $\mathbf{u}_j$. 

Specifically, $\sigma_{E}^{2}(\mathbf{u})$ is defined by:

$$
\sigma_E^2(\mathbf{u}) = 
    - \sum_{i=1}^n \sum_{j=1}^n 
\lambda_i \lambda_j 
\gamma(\mathbf{u}_i - \mathbf{u}_j) +
    2 \sum_{j=1}^n \lambda_i \gamma(\mathbf{u}_i - \mathbf{u}),
\tag{11a}
\label{eq:kriging_sigma}
$$

which can be expressed in matrix notation as it relates to $\eqref{eq:kriging_matrix}$:

$$
\sigma_E^2(\mathbf{u}) = 
-\lambda^\text{T} \mathbf{\Gamma}_{n \times n} + 
2\lambda^\text{T}\gamma_\mathbf{u},
\tag{11b}
$$

where:

$$
\begin{equation}
\lambda =
\begin{bmatrix}
    \lambda_1 \\
    \lambda_2 \\
    \vdots \\
    \lambda_n
\end{bmatrix},
\quad
\gamma_\mathbf{u} =
\begin{bmatrix}
    \gamma(\mathbf{u}_1 - \mathbf{u}) \\
    \gamma(\mathbf{u}_2 - \mathbf{u}) \\
    \vdots \\
    \gamma(\mathbf{u}_n - \mathbf{u})
\end{bmatrix}.
\tag{11c}
\end{equation}
$$

### Solving for $\lambda$

The constrained minimization problem, minimizing $\sigma_E^2$ subject to $\sum\lambda_b = 1$, is solved using the **Lagrange Multiplier technique**. This is accomplished by introducing a the Lagrange multiplier, $\mu$, which enforces the constraint during the minimization. The objective function $\eqref{eq:kriging_sigma}$ is combined with the constraint using the Lagrangian function:

$$
\mathcal{L}(\lambda_1, \ldots, \lambda_n, \mu) = \sigma_E^2(\mathbf{u}) + 2\mu \left(1 - \sum_{j=1}^n \lambda_j\right).
\tag{12}
$$

To find the minimum of this Lagrangian, partial derivatives are taken with respect to each weight $\lambda_b$ and with respect to the Lagrange multiplier $\mu$. Setting these derivatives to zero yields:

$$ 
\frac{\partial \mathcal{L}}{\partial \lambda_b} = 2 \sum_{j=1}^n \lambda_j \gamma(\mathbf{u}_j - \mathbf{u}_b) - 2 \gamma(\mathbf{u}_b - \mathbf{u}) - 2\mu = 0,
\tag{13a} 
$$

which simplifies to:

$$
\sum_{j=1}^n \lambda_j \gamma(\mathbf{u}_j - \mathbf{u}_b) + \mu = \gamma(\mathbf{u}_b - \mathbf{u}).
\label{eq:partialdev_simple}
\tag{13b}
$$

The partial derivative with respect to $\mu$ returns the constraint:

$$
\frac{\partial \mathcal{L}}{\partial \mu} = 2 \left(1 - \sum_{j=1}^n \lambda_j\right) = 0,
\label{eq:constraint}
\tag{13c}
$$

which is $\sum\limits_{j=1}^n \lambda_j = 1$. Using both $\eqref{eq:partialdev_simple}$ and $\eqref{eq:constraint}$ forms the system of equations:

$$
\begin{equation}
\begin{cases}
    \sum\limits_{j=1}^n \lambda_j(\mathbf{u})\gamma(\mathbf{u}_j - \mathbf{u}_i) + \mu = \gamma(\mathbf{u}_i - \mathbf{u}) & \text{for } i = 1, \ldots, n \\
    \sum\limits_{j=1}^n \lambda_j(\mathbf{u}) = 1
\end{cases}.
\tag{14}
\label{eq:lagrange_sys}
\end{equation}
$$

This system can be written in matrix form as:

$$
\mathbf{\Gamma} \hat{\mathbf{\lambda}} = \hat{\mathbf{\gamma}}_\mathbf{u},
\tag{15a}
$$

where:

$$
\begin{equation}
\underbrace{
    \begin{bmatrix}
        \gamma_{1,1} & \cdots & \gamma_{1,n} & 1 \\
        \vdots & \ddots & \vdots & \vdots \\
        \gamma_{n,1} & \cdots & \gamma_{n,n} & 1 \\
        1 & \cdots & 1 & 0
    \end{bmatrix}
}_{\mathbf{\Gamma}}
\underbrace{
    \begin{bmatrix}
        \lambda_1 \\
        \vdots \\
        \lambda_n \\
        \mu
    \end{bmatrix}
}_{\hat{\mathbf{\lambda}}} =
\underbrace{
    \begin{bmatrix}
        \gamma(\mathbf{u}_1 - \mathbf{u}) \\
        \vdots \\
        \gamma(\mathbf{u}_n - \mathbf{u}) \\
        1
    \end{bmatrix}
}_{\mathbf{\hat{\gamma}}_\mathbf{u}}.
\tag{15b}
\end{equation}
$$

In these equations, $\gamma_{i,j} = \gamma(\mathbf{u}_i - \mathbf{u}_j)$ is the semivariance between data points $i$ and $j$, and $\gamma_{i,\mathbf{u}} = \gamma(\mathbf{u}_i - \mathbf{u})$ is the semivariance between data point $i$ and the target location $\mathbf{u}$. The final row, final column, and the bottom-right zero in $\mathbf{\Gamma}$ are generated by $\eqref{eq:lagrange_sys}$  to enforce the unbiasedness constraint $\sum \lambda_i = 1$.

Solving this system yields the kriging weights $\lambda_1, \ldots, \lambda_n$ and $\mu$. The Lagrange multiplier is essential for enforcing the unbiasedness constraint and appears as an additional unknown in the system. This approach ensures that the kriging predictor is both unbiased and has minimum variance among all linear unbiased estimators. In practical computation, the system $\mathbf{\Gamma} \hat{\mathbf{\lambda}} = \hat{\mathbf{\gamma}}_\mathbf{u}$ is solved for the vector of unknowns, which includes both $\hat{\lambda}$ and $\mu$. The solution provides the optimal $\hat{\lambda}$ to be applied to $z(\mathbf{u}_b)$ , as well as the value of the $\mu$ that enforces the unbiasedness constraint. 

## Kriging prediction and uncertainty

Once the system for the kriging weights $\lambda_1, \ldots, \lambda_n$ and $\mu$ has been solved, the actual interpolation at the target location $\mathbf{u}$, $\mathbf{z}^*(\mathbf{u})$, is performed using the predictor equation defined previously $\eqref{eq:kriging_pred}$. Each prediction location will generally yield a different set of weights, reflecting the spatial configuration and values of the observed data. The kriging variance, or mean squared prediction error, quantifies the expected uncertainty of the prediction at $\mathbf{u}$. As shown above, it is computed as:

$$
\sigma_K^2(\mathbf{u}) = \sum_{i=1}^n \lambda_i(\mathbf{u})\, \gamma(\|\mathbf{u}_i - \mathbf{u}\|) + \mu - \gamma(0),
\tag{16}
$$

where $\gamma(\|\mathbf{u}_i - \mathbf{u}\|)$ is the semivariogram value between each data location and the prediction location, $\mu$ is the Lagrange multiplier from the kriging system, and $\gamma(0)$ is the process variance. This variance provides a local measure of the reliability of the prediction: lower values indicate higher confidence, while higher values suggest greater uncertainty, often due to sparse or distant data.

In addition to the kriging variance, it is often useful to report the sample variance or the coefficient of variation (CV) for each prediction. The coefficient of variation is defined as

$$
\text{CV}(\mathbf{u}) = \frac{\sqrt{\sigma_K^2(\mathbf{u})}}{|\mathbf{z}^*(\mathbf{u})|}.
\tag{17}
$$

The CV provides a dimensionless measure of relative uncertainty, allowing for comparison across locations with different magnitudes of predicted values.

When kriging predictions are made over a spatial mesh or grid, it is often necessary to aggregate results to obtain area-weighted means, totals, or global uncertainty measures. The area-weighted mean of predictions over $M$ grid cells with areas $A_m$ is:

$$
\bar{\mathbf{z}}_{\text{area}}^* = \frac{\sum\limits_{m=1}^M A_m z^*(\mathbf{u}_m)}{\sum\limits_{m=1}^M A_m}.
\tag{18}
$$

The area-weighted variance can be computed as:

$$
\mathbb{V}_{\text{area}} = \frac{\sum\limits_{m=1}^M A_m^2 \sigma_K^2(\mathbf{u}_m)}{\left(\sum\limits_{m=1}^M A_m\right)^2}.
\tag{19}
$$

These aggregated quantities are useful for summarizing spatial predictions and uncertainties over regions of interest, such as survey strata or management zones.


In summary, the ordinary kriging workflow consists of: (1) modeling the spatial dependence via the semivariogram, (2) solving the kriging system for each prediction location to obtain the weights and Lagrange multiplier, (3) computing the interpolated value as a weighted sum of observed data, and (4) quantifying the prediction uncertainty using the kriging variance. This approach ensures that predictions are both unbiased and have minimum variance among all linear estimators, provided the model assumptions are satisfied.

## Practical considerations

- Variogram/covariance choice: The model must capture short-range behaviour (near the origin) and an appropriate range scale; mismatches bias weights and uncertainty estimates.[^cressie]  
- Conditioning: The augmented matrix in can be ill-conditioned in some sampling geometries; numerical stabilization or mild regularization is often necessary in practice.  
- Neighborhood definition: Limiting the set of samples used for a given prediction trades off locality against numerical stability and computational cost; this is a modeling choice rather than a theoretical requirement.  
- Extrapolation: Predictions made outside the effective range of spatial dependence should be treated with caution; reported kriging variances may understate true uncertainty when extrapolation is extensive.[^rivoirard]

## References

[^cressie]: Cressie, N. (1993). *Statistics for Spatial Data*. Wiley.  
[^chiles]: Chilès, J.-P. & Delfiner, P. (2012). *Geostatistics: Modeling Spatial Uncertainty*. Wiley.  
[^journel]: Journel, A.G. & Huijbregts, C.J. (1978). *Mining Geostatistics*. Academic Press.  
[^rivoirard]: Rivoirard, J. et al. (2000). *Geostatistics for Estimating Fish Abundance*. Blackwell Science.