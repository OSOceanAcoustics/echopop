(variogram-theory)=
# Characterizing spatial correlation

## Stationarity and the semivariogram

In geostatistics, a **spatial process** or **random field** is a collection of random variables:

<a id="eq:randomz"></a>

$$
\{ Z(\mathbf{x}) : \mathbf{x} \in D \subset \mathbb{R}^d\},
\tag{1}
\label{eq:randomz}
$$

where $D \subset \mathbb{R}^d$ is the spatial domain, $d \in \mathbb{N}$ is the dimension of the space, and $Z(\mathbf{x})$ is a random variable that represents the value of some spatially varying quantity (e.g. temperature) at location $\mathbf{x}$. Importantly, the term *process* refers to this entire stochastic system, not a single observation at any specific location. A spatial process is said to be **second-order stationary** if its expected value is constant and its covariance depends only on the separation between two locations:

$$
\mathbb{E}[Z(\mathbf{x})] = m, \tag{2a}
$$

$$
\mathbb{E}\!\left[(Z(\mathbf{x}) - m)(Z(\mathbf{x}+\mathbf{h}) - m)\right] = 
\text{Cov}\big(Z(\mathbf{x}),Z(\mathbf{x}+\mathbf{h})\big) =
\mathbb{C}(h), 
\tag{2b}
$$

where $m$ is a constant expected process mean and $\mathbb{C}(\cdot)$ the covariance function. Here $\mathbf{h}$ denotes the lag vector (direction and magnitude) between two locations and $h$ denotes the Euclidean length (scalar lag distance). Either $\mathbf{h}$ or $h$ can be used depending on whether directionality matters, which is discussed later ([**isotropy and anisotropy**](#isotropy-and-anisotropy)). For a specific pair of separate locations indexed by $i$ and $j$ ($i, j \in \{1, \dots, n_s \}$, where $n_s$ is the number of samples) define $\mathbf{h}$ and $h$ by:

$$
\mathbf{h} = \mathbf{x}_j - \mathbf{x}_i, 
\qquad
h = \| \mathbf{h} \| = \| \mathbf{x}_j - \mathbf{x}_i \| = 
\sqrt{\sum_{k=1}^d \big(x_{j,k} - x_{i,k}\big)^2},
\tag{3}
$$

where $x_{i,k}$ denotes the $k$‑th coordinate of $\mathbf{x}_i$ and $k = 1 \dots d$ indexes spatial coordinates (e.g., $d=1$ for a line, $d=2$ for a plane and $d=3$ for three‑dimensional space). The specific pairwise scalar distance used in empirical sums and binning is written $h_{ij}=\|\mathbf{x}_j-\mathbf{x}_i\|$, whereas plain $h$ denotes a generic lag argument in theoretical expressions such as $\gamma(h)$ or $\mathbb{C}(h)$.

:::{admonition} Theoretical vs. empirical expressions
:class: note
While $Z(\mathbf{x})$ denotes the random field, the quantity $z(\mathbf{x})$ can be used to express a realized observation at location $\mathbf{x}$. This distinction makes clear whether an expression is a theoretical statement or an empirical computation.
:::

This information helps characterize spatial dependence by using the **semivariogram**, which describes how the average squared difference between values changes with distance. The **semivariance** at lag $h$ is defined as:

$$
\gamma(h) = \mathbb{C}(0) - \mathbb{C}(h),
\tag{4}
\label{eq:base_semivariance}
$$

where $\mathbb{C}(0)$ is the variance at zero lag (the process variance). In this context, the semivariogram and covariance function contain equivalent information about spatial correlation. In this document, the scalar $h$ denotes a generic lag distance.

## Isotropy and anisotropy

A spatial process is said to be **isotropic** if its spatial dependence structure is the same in all directions; that is, the covariance and semivariogram functions depend only on the distance between locations, not on the direction. Mathematically, this means that the semivariogram and covariance are functions of $h$, hence $\eqref{eq:base_semivariance}$. This assumption simplifies variogram analysis, as all pairs of points separated by the same distance are treated equivalently, regardless of their orientation.

In contrast, a process is **anisotropic** if its spatial dependence varies with direction. In this case, the covariance and semivariogram depend on both the magnitude and orientation of $\mathbf{h}$:

$$
\gamma(\mathbf{h}) = \mathbb{C}(\mathbf{0}) - \mathbb{C}(\mathbf{h}).
\tag{5}
$$

Directional (azimuthal) variogram analysis is required to detect and model anisotropy. This involves grouping pairs of locations not only by lag distance but also by their relative orientation (azimuth angle), allowing the spatial correlation structure to be characterized as a function of both distance and direction.

## The intrinsic model

However, real-world datasets often do not fulfill the requirements for full second-order stationarity. The **intrinsic hypothesis** provides a weaker condition by requiring only that the *increments* of $Z(\mathbf{x})$, $Z(\mathbf{x}) - Z(\mathbf{x}+\mathbf{h})$, be second-order stationary. Under this weaker condition, the semivariogram is now defined as:

$$
\gamma(h) = \frac{1}{2} \mathbb{E}\left[ (Z(\mathbf{x}) - Z(\mathbf{x}+\mathbf{h}))^2 \right],
\tag{6}
$$

which means that, for a given $\mathbf{h}$, all possible pairs of locations in the spatial domain separated by $\mathbf{h}$ are considered. The expectation $\mathbb{E}[\cdot]$ then represents the average of these squared differences over all such pairs. Thus, the intrinsic model generalizes the covariance-based approach, allowing spatial correlation to be characterized through the semivariance even when the covariance function is not well-defined. In practice, however, $Z(\mathbf{x})$ is not known everywhere; only a finite set of observations at sampled locations is available. 

To estimate the semivariogram from data, the theoretical expectation is replaced with an average over all observed pairs of points separated by lag $h$. This leads to the **empirical semivariogram** estimator:

$$
\hat{\gamma}(h) = \frac{1}{2N(h)} \sum_{i<j:\; h_{ij}\approx h} \bigl(z(\mathbf{x}_j)-z(\mathbf{x}_i)\bigr)^2,
\tag{7}
$$

where $z(\mathbf{x}_i)$ and $z(\mathbf{x}_j)$ denote the observed (realized) values of the variable at sample locations $\mathbf{x}_i$ and $\mathbf{x}_j$, respectively. Each $\mathbf{x}_i$ (or $\mathbf{x}_j$) is a location in the sampling set; $z(\mathbf{x}_i)$ is the concrete measurement recorded at that location (scalar unless a multivariate observation is specified). This estimator provides a practical method for computing the semivariogram from data, enabling empirical characterization of spatial dependence.

### Practical notes on empirical computation

When computing the empirical semivariogram, the choice of lag bins is important. Lag bins are typically defined by dividing the range of observed pairwise distances into intervals of fixed width, or by using quantiles to ensure a similar number of pairs per bin. The bin width and maximum lag should be chosen to balance resolution and statistical reliability; bins with too few pairs can yield unstable estimates. In practice, a tolerance is used so that all pairs with $h_{ij}$ within a specified interval around $h$ are included in the bin for $h$. Uneven sampling or missing data may result in some bins having very few pairs, which can be mitigated by adjusting bin widths or excluding sparse bins from model fitting. Formally, $N(h)$ is the count of unique pairs falling in the bin for $h$:

$$
N(h)=\sum_{i<j} \mathbf{1}\{\,h_{ij}\in [h-\Delta h/2,\;h+\Delta h/2]\,\},
\tag{8}
$$

where $\Delta h$ is the bin width and $\mathbf{1}\{\cdot\}$ is the indicator function.


## Standardizing empirical variograms

The classic empirical semivariogram assumes homoscedasticity (homogeneous variance) across space. Many real datasets exhibit heteroscedasticity, where local variance changes with  location due to patchiness, measurement error, or sampling design. This can bias or destabilize semivariogram estimates. To reduce sensitivity to local variance differences and outliers, the empirical semivariogram from $\eqref{eq:empirical_semivariogram}$ is normalized by the product of their local standard deviations:

$$
\hat{\gamma}_{\text{std}}(h) = \frac{1}{2N(h)} \sum_{i<j:\; h_{ij}\approx h} \frac{[z(\mathbf{x}_j) - z(\mathbf{x}_i)]^2}{\sigma(\mathbf{x}_j) \cdot \sigma(\mathbf{x}_i)},
\tag{9}
$$

where $\sigma(\mathbf{x}_i)$ and $\sigma(\mathbf{x}_j)$ are local standard deviations at $\mathbf{x}_i$ and $\mathbf{x}_j$. This standardization improves robustness to local variance differences and outliers, providing more reliable estimates of spatial correlation in heterogeneous datasets.

## Semivariogram features and interpretation

Theoretical semivariogram models are parametric functions designed to capture the essential features of spatial correlation observed in empirical semivariograms. These models are characterized by several key parameters:

- **Nugget ($C_0$)**: The value of the semivariogram at zero lag. The nugget represents measurement error or micro-scale variation that occurs at distances smaller than the sampling interval. A large nugget indicates substantial noise or unresolved spatial structure.

- **Sill ($C$)**: The plateau value that the semivariogram approaches at large lag distances. The sill reflects the total variance of the process. In bounded models, the semivariogram levels off at the sill; in unbounded models, it may increase indefinitely.

- **Range ($a$)**: The lag distance at which the semivariogram reaches the sill (or a specified fraction of the sill, such as 95%). The range indicates the scale of spatial correlation: pairs of points separated by less than the range are correlated, while those farther apart are effectively uncorrelated.

- **Model shape**: The functional form of the model determines how quickly spatial correlation decays with distance. For example, some models rise linearly near the origin and flatten at the sill (spherical), others rise rapidly and approach the sill asymptotically (exponential), and some are very smooth near the origin (Gaussian). The choice of model shape reflects the expected smoothness and continuity of the underlying spatial process.

- **Anisotropy**: In some cases, spatial correlation may vary with direction. Anisotropic models allow the range or sill to depend on the orientation of the lag vector, capturing directional effects in the data.

Interpreting these features provides insight into the spatial structure of the data, the scale of spatial dependence, and the presence of measurement error or micro-scale variability.

A comprehensive list of theoretical variogram models, including their mathematical equations and descriptive properties, is provided in the **[semivariogram equations](semivariogram_eq.md)** reference document. Readers should consult that document for details on specific models, their parameterizations, and guidance on model selection for different types of spatial processes.

### Diagnostics and interpretation

After fitting a theoretical variogram model, it is important to assess the quality of the fit. Visual inspection of the empirical and fitted curves can reveal systematic deviations, such as underfitting at short or long lags, or failure to capture anisotropy. Quantitative diagnostics, such as cross-validation or analysis of kriging residuals, can help evaluate whether the chosen model adequately represents spatial dependence. Poor fits may indicate the need for alternative model shapes, inclusion of anisotropy, or further data cleaning. Interpreting the fitted parameters in the context of the data and sampling design is essential for drawing reliable conclusions about spatial structure.

## Semivariogram model fitting

After estimating the empirical semivariogram from spatial data, a theoretical variogram model is fit to the empirical values to obtain a smooth, interpretable representation of spatial correlation. Theoretical models are parametric functions of the form:

$$
\gamma(h) = C_0 + C_1 \cdot \mathscr{f}(h; \theta),
\tag{10}
$$

where $C_0$ is the nugget (the value at zero lag), $C_1$ is the partial sill ($C - C_0$), and $\mathscr{f}(h; \theta)$ is a correlation function that depends on lag $h$ and model-specific parameters $\theta$ (e.g., range, smoothness, periodicity). The process of fitting a variogram model involves finding the set of parameters $\theta$ that best match the empirical semivariogram. This is typically achieved by minimizing the weighted sum of squared differences between the empirical semivariogram values $\gamma_\text{empirical}(h_i)$ and the theoretical model $\gamma_\text{model}(h_i; \theta)$ at each lag $h_i$:

$$ 
\min_{\theta} \sum_{b=1}^{n} w_b \left[ \gamma_\text{empirical}(h_b) - \gamma_\text{model}(h_b; \theta) \right]^2,
\tag{11}
$$

where $n$ is the number of lag bins and $b$ indexes those lag bins ($b=1, \dots, n$). The weights $w_b$ are associated with each lag bins. Values of $w_b$ are often chosen to be proportional to the number of data pairs in the $b$‑th bin, $N(h_b)$, so that lag bins with more data have greater influence on the fit. When no standardization is applied, these weights are defined as:

$$
w_b = N(h_b).
\tag{12}
$$

In contrast, when $w_b$ are normalized such that they sum to 1:
 
$$
\hat{w}_b = \frac{N(h_b)}{\sum\limits_{b=1}^{n} N(h_b)}.
\tag{13}
$$

## Assumptions and limitations

Variogram analysis assumes at least intrinsic stationarity of increments and, for second‑order stationarity, a constant mean. If the mean exhibits spatial trend, detrend or model the mean prior to variogram estimation. Variogram estimation is sensitive to outliers, edge effects, and irregular sampling; handle these issues explicitly and document the choices made. Standardized estimators, robust estimators (e.g. Cressie–Hawkins), and careful binning can mitigate some practical problems, but they do not replace the need to evaluate assumptions and diagnostics.


## References and further reading

For more details on variogram theory, estimation, and applications, see:

- Cressie, N. (1993). *Statistics for Spatial Data*. Wiley.
- Journel, A. G., & Huijbregts, C. J. (1978). *Mining Geostatistics*. Academic Press.
- Webster, R., & Oliver, M. A. (2007). *Geostatistics for Environmental Scientists*. Wiley.

These texts provide comprehensive coverage of spatial statistics, variogram modeling, and practical guidance for geostatistical analysis.