(variogram-theory)=
# Characterizing spatial correlation

## Stationarity and the semivariogram

In geostatistics, a **spatial process** or **random field** is a collection of random variables,

$$
\{ Z(\mathbf{x}) : \mathbf{x} \in D \subset \mathbb{R}^d\},
\tag{1}
$$

where $D \subset \mathbb{R}^d$ is the spatial domain, $d \in \mathbb{N}$ is the dimension of the space, and $Z(\mathbf{x})$ is a random variable that represents the value of some spatially varying quantity (e.g. temperature) at location $\mathbf{x}$. Importantly, the term *process* refers to this entire stochastic system, not a single observation at any specific location. A spatial process is said to be **second-order stationary** if its expected value is constant and its covariance depends only on the separation between two locations:

$$
\mathbb{E}[Z(\mathbf{x})] = m, \tag{2a}
$$

$$
\mathbb{E}\!\left[(Z(\mathbf{x}) - m)(Z(\mathbf{x}+\mathbf{h}) - m)\right] = C(h), \tag{2b}
$$

where $m$ is a constant expected process mean, $\mathbf{h}$ is the lag vector separating two locations, $h$ is the lag distance, and $C(h)$ is the covariance function that describes how values at two locations, ($x_i, y_i$) and ($x_j, y_j$), are correlated as a function of distance. The Euclidean norm of $\mathbf{h}$ defines $h$ between $\mathbf{x}_i$ and $\mathbf{x}_j$ ($h_{ij}$):

$$
h_{ij} = |\mathbf{x}_i - \mathbf{x}_j| = \sqrt{(x_i - x_j)^2 + (y_i - y_j)^2}.
\tag{3}
$$

This information helps characterize spatial dependence by using the **semivariogram**, which describes how the average squared difference between values changes with distance. The **semivariance** at lag $h$ is defined as:

$$
\gamma(h) = C(0) - C(h),
\tag{4}
\label{eq:base_semivariance}
$$

where $C(0)$ is the variance at zero lag (the process variance). In this case, the semivariogram and covariance function contain equivalent information about spatial correlation.

## Isotropy and anisotropy

A spatial process is said to be **isotropic** if its spatial dependence structure is the same in all directions; that is, the covariance and semivariogram functions depend only on the distance between locations, not on the direction. Mathematically, under isotropy, the semivariogram and covariance are functions of the scalar lag $h = \|\mathbf{h}\|$ via $\eqref{eq:base_semivariance}$. This assumption simplifies variogram analysis, as all pairs of points separated by the same distance are treated equivalently, regardless of their orientation.

In contrast, a process is **anisotropic** if its spatial dependence varies with direction. In this case, the covariance and semivariogram depend on both the magnitude and orientation of the lag vector $\mathbf{h}$:

$$
\gamma(\mathbf{h}) = C(\mathbf{0}) - C(\mathbf{h}).
\tag{5}
$$

Directional (azimuthal) variogram analysis is required to detect and model anisotropy. This involves grouping pairs of locations not only by lag distance but also by their relative orientation (azimuth angle), allowing the spatial correlation structure to be characterized as a function of both distance and direction.

## The intrinsic model

However, real-world datasets often do not fulfill the requirements for full second-order stationarity. The **intrinsic hypothesis** provides a weaker condition by requiring only that the *increments* of $Z(\mathbf{x})$, $Z(\mathbf{x}) - Z(\mathbf{x}+\mathbf{h})$, be second-order stationary. Under this weaker condition, the semivariogram is now defined as:

$$
\gamma(h) = \frac{1}{2} \mathbb{E}\left[ (Z(\mathbf{x}) - Z(\mathbf{x}+\mathbf{h}))^2 \right],
\tag{5}
$$

which means that, for a given $\mathbf{h}$, all possible pairs of locations in the spatial domain separated by $\mathbf{h}$ are considered. The expectation $\mathbb{E}[\cdot]$ then represents the average of these squared differences over all such pairs. Thus, the intrinsic model generalizes the covariance-based approach, allowing spatial correlation to be characterized through the semivariance even when the covariance function is not well-defined. In practice, however, $Z(\mathbf{x})$ is not known everywhere; only a finite set of observats at sampled locations is available. 

To estimate the semivariogram from data, the theoretical expefctation is replaced with an average over all observed pairs of points separated by lag $h$. This leads to the **empirical semivariogram** estimator:

$$
\hat{\gamma}(h) = \frac{1}{N(h)} \sum_{i=1}^{N(h)} \frac{1}{2} \left(z(s_i) - z(s_i+h)\right)^2,
\tag{6a}
\label{eq:empirical_semivariogram}
$$

where $N(h)$ is the number of pairs of observed locations separated by $h$, and $z(s_i)$ and $Z(s_i+h)$ are the observed values at locations $s_i$ and $s_i+h$. This estimator provides a practical method for computing the semivariogram from data, enabling empirical characterization of spatial dependence. Alternatively, this can expressed in terms of $\mathbf{x}_i$ and $\mathbf{x}_j$:

$$ 
\hat{\gamma}(h) = \frac{1}{2N(h)} \sum_{i,j: h_{ij} \approx h} \left(z(\mathbf{x}_i) - z(\mathbf{x}_j)\right)^2.
\tag{6b}
$$

## Standardizing empirical variograms

The classic empirical semivariogram assumes homoscedasticity (homogenous variance) across space. Many real-world datasets exhibit heteroscedatisticity, where local variance changes with locations due to patchiness, measurement error, or sampling design. This can bias or destabilize semivariogram estimates. To address this, the empirical semivariogram from $\eqref{eq:empirical_semivariogram}$ normalizes the squared differences between paired observations by the product of their local standard deviations:

$$
\hat{\gamma}_{\text{std}}(h) = \frac{1}{2N(h)} \sum_{i,j: h_{ij}} \frac{[Z(\mathbf{x}_i) - Z(\mathbf{x}_j)]^2}{\sigma(\mathbf{x}_i) \cdot \sigma(\mathbf{x}_j)},
\tag{7}
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


## Semivariogram model fitting

After estimating the empirical semivariogram from spatial data, a theoretical variogram model is fit to the empirical values to obtain a smooth, intrpretable representation of spatial correlation. Theoretical models are parametric functions of the form:

$$
γ(h) = C_0 + C_1 \cdot f(h; θ),
\tag{8}
$$

where $C_0$ is the nugget (the value at zero lag), $C_1$ is the partial sill ($C - C_0$), and $f(h; \theta)$ is a correlation function that depends on lag $h$ and model-specific parameters $\theta$ (e.g. range, smoothness, periodicity). The process of fitting a variogram model involves finding the set of parameters $\theta$ that best match the empirical semivariogram. This is typically achieved by minimizing the weighted sum of squared differences between the empirical semivariogram values $\gamma_\text{empirical}(h_i)$ and the theoretical model $\gamma_\text{model}(h_i; \theta)$ at each lag $h_i$:

$$ 
\min_{\theta} \sum_{i=1}^{n} w_i \left[ \gamma_\text{empirical}(h_i) - \gamma_\text{model}(h_i; \theta) \right]^2,
\tag{9}
$$

where $n$ is the number of lag bins and $w_i$ are weights. Values of $w_i$ areoften chosen to be proportional to the number of data pairs at each lag $h_i$, so that lag bins with more data have greater influence on the fit. When no standardization is applied, these weights are defined as:

$$
w_i = N(h_i).
\tag{10}
$$

In contrast, when $w_i$ are normalized such that they sum to 1:

$$
\hat{w}_i = \frac{N(h_i)}{\sum\limits_{j=1}^{n} N(h_j)}.
\tag{11}
$$
