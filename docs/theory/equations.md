---
orphan:
---

# Equations for internal referencing

(eq-16)=
:::{admonition} TS-length relationship
$$
    \textit{TS} = m\log_{10}L+b
$$

See {ref}`Equation 1.6 <intext_eq_16_md>` for more details.

:::

(eq-110)=
:::{admonition} Relationship between volumetric backscatter and animal number density (animals m<sup>-3</sup>)
$$
    s_\text{v} = \rho_\text{v} \left< \sigma_\text{bs} \right>
$$

See {ref}`Equation 1.10 <intext_eq_110_md>` for more details.

:::

(eq-111)=
:::{admonition} Vertical integration of volumetric backscatter in the water column
$$
    s_\text{a} = 
        \int\limits_{z_1}^{z_2} s_\text{v} dz =
            s_\text{v} H
$$

See {ref}`Equation 1.11 <intext_eq_111_md>` for more details.

:::

(eq-113)=
:::{admonition} Vertical integration of volumetric backscatter in the water column
$$
    \rho_\text{a}(x) = \frac{s_\text{a}(x)}{\left< \sigma_\text{bs} \right>}
$$

See {ref}`Equation 1.13 <intext_eq_113_md>` for more details.

:::

(eq-115)=
:::{admonition} Areal number density (animals nmi<sup>-2</sup>)
$$
    \rho_\text{A} =
        \frac{s_\text{A}}{4 \pi \left< \sigma_\text{bs} \right>} =
            \frac{s_\text{A}}{\sigma_\text{sp}}
$$

See {ref}`Equation 1.15 <intext_eq_115_md>` for more details.

:::

(eq-21b)=
:::{admonition} Average animal weights
$$
    \left< w \right> = \frac{\sum\limits_{j=1}^N w_j}{N}
$$

See {ref}`Equation 2.1b <intext_eq_21_md>` for more details.

:::

(eq-23)=
:::{admonition} Normalized length frequency distribution
\begin{align*}
    \tilde{\mathbf{L}} =
    \left[
        \begin{split}
            \tilde{L}&_1 \\ \tilde{L}&_2 \\ \tilde{L}&_3 \\ \vdots&
        \end{split}
    \right].
\end{align*}

See {ref}`Equation 2.3 <intext_eq_23_md>` for more details.

:::

(eq-25)=
:::{admonition} Weight-length distribution
\begin{align*}
    \mathbf{w} =
    \left[
        \begin{split}
            w&_1 \\ w&_2 \\ w&_3 \\ \vdots&
        \end{split}
    \right].
\end{align*}

See {ref}`Equation 2.5 <intext_eq_25_md>` for more details.

:::

(eq-26)=
:::{admonition} Weight-length log-linear relationship
$$
    \log_{10}(\hat{w}) =
        \log_{10}(\hat{a}) + \hat{b} \log_{10}(L)
$$

See {ref}`Equation 2.6 <intext_eq_26_md>` for more details.

:::

(eq-27)=
:::{admonition} Weight estimate for length bin $\ell$
$$
    w_\ell =
    \left[
        10^{\hat{a}} {(L_{\ell}^{*})}^{\hat{b}}
    \right]
    L_\ell
$$

See {ref}`Equation 2.6 <intext_eq_26_md>` for more details.

:::

(eq-29)=
:::{admonition} Areal biomass density (kg nmi<sup>-2</sup>) for each sex
$$
    \rho_{B; s} = \rho^i_{B; s}(x,y) = \rho_A(x,y) (\tilde{\mathbf{L}}^i_s)^\top \mathbf{w}^i_s
$$

See {ref}`Equation 2.9 <intext_eq_29_md>` for more details.

:::

(eq-210)=
:::{admonition} Sex-specific normalized length frequency distribution
\begin{align*}
    \tilde{\mathbf{L}}^i_s &=
    \left[
        \begin{split}
            \tilde{L}&^i_{s,1} \\ \tilde{L}&^i_{s,2} \\ \tilde{L}&^i_{s,3} \\ &\vdots
        \end{split}
    \right],
    \\[2ex]
        \sum_{s,\ell} \tilde{L}^i_{s,\ell} &= 1.
\end{align*}

See {ref}`Equation 2.10 <intext_eq_210_md>` for more details.

:::

(eq-211)=
:::{admonition} Normalized weight-length distribution for each sex
\begin{align*}
    \tilde{\mathbf{w}}^i_s =
        \left[
        \begin{split}
            \tilde{w}&^i_{s,1} \\ \tilde{w}&^i_{s,2} \\ \tilde{w}&^i_{s,3} \\ &\vdots
        \end{split}
    \right].
\end{align*}

See {ref}`Equation 2.11 <intext_eq_211_md>` for more details.

:::

(eq-212)=
:::{admonition} Normalized weight-length distribution bins
\begin{align*}
        \tilde{w}^i_{s,\ell} &= 
            \frac{\mathbf{w}^i_s}{\sum\limits_{\ell} w^i_{s,\ell}} =
                \frac{w^i_\ell}{\sum\limits_{s,\ell} w^i_{s,\ell}},
    \\[2ex]
        \sum\limits_{s,\ell} \tilde{w}^i_{s,\ell} &= 1. 
\end{align*}

See {ref}`Equation 2.12 <intext_eq_212_md>` for more details.

:::

(eq-215)=
:::{admonition} Transect-specific mean density
$$
    \hat{z}^{\,t} = 
        \frac{\sum\limits_{k \in t} z(x^k)}{d^t}
$$

See {ref}`Equation 2.15 <intext_eq_215_md>` for more details.

:::

(eq-216)=
:::{admonition} Distance-based transect weights
$$
    \tau^t = 
        \frac{d^t}{\frac{1}{n^i} \sum\limits_{t \in i} d^t}
$$

See {ref}`Equation 2.16 <intext_eq_216_md>` for more details.

:::

(eq-217)=
:::{admonition} Stratum-specific mean density

$$
    \hat{z}^{i} = \frac{1}{n^i} \sum\limits_{t \in i} \tau^t \hat{z}^{t}
$$

See {ref}`Equation 2.17 <intext_eq_217_md>` for more details.

:::

(eq-218)=
:::{admonition} Survey mean density

$$
    \hat{z} = \frac{\sum\limits_{i} A_i \hat{z}^i}{\sum\limits_i A_i}
$$

See {ref}`Equation 2.18 <intext_eq_218_md>` for more details.

:::

(eq-219)=
:::{admonition} Within-stratum survey variance

$$
    \mathbb{V}(\hat{z}^i) =
            \frac{\sum\limits_{t \in i} \tau^t (\hat{z}^t - \hat{z}^i)^2}{n^i(n^i - 1)}
$$

See {ref}`Equation 2.19 <intext_eq_219_md>` for more details.

:::

(eq-220)=
:::{admonition} Area-weighted survey variance

$$
    \mathbb{V}(\hat{z}) =
            \frac{\sum\limits_i (A^i)^2\, \mathbb{V}(\hat{z}^i)}{\left( \sum\limits_i A^i \right)^2}
$$

See {ref}`Equation 2.20 <intext_eq_220_md>` for more details.

:::

(eq-31)=
:::{admonition} Random field (or spatial process)

$$
    \{ Z(\mathbf{x}) : \mathbf{x} \in D \subset \mathbb{R}^d\}
$$

See {ref}`Equation 3.1 <intext_eq_31_md>` for more details.

:::

(eq-34)=
:::{admonition} Semivariance

$$
    \gamma(h) = \mathbb{C}(0) - \mathbb{C}(h)
$$

See {ref}`Equation 3.4 <intext_eq_34_md>` for more details.

:::

(eq-37)=
:::{admonition} Empirical semivariogram

$$
    \hat{\gamma}(h) = \frac{1}{2N(h)} \sum_{i<j:\; h_{ij}\approx h} \bigl(z(\mathbf{x}_j)-z(\mathbf{x}_i)\bigr)^2
$$

See {ref}`Equation 3.7 <intext_eq_37_md>` for more details.

:::

(eq-319)=
:::{admonition} Kriging linear predictor

$$
    \mathbf{z}^*(\mathbf{u}) = 
    \sum_{b=1}^n \lambda_b(\mathbf{u}) z(\mathbf{u}_b) + 
        \left[ 
            1 - \sum_{b=1}^n \lambda_b(\mathbf{u}) 
        \right] \mathbf{m}
$$

See {ref}`Equation 3.19 <intext_eq_319_md>` for more details.

:::

(eq-321)=
:::{admonition} Ordinary Kriging linear predictor

$$
    \mathbf{z}^*(\mathbf{u}) = 
    \sum_{b=1}^n \lambda_b(\mathbf{u}) z(\mathbf{u}_b) + 
        \left[ 
            1 - \sum_{b=1}^n \lambda_b(\mathbf{u}) 
        \right] \mathbf{m}
$$

See {ref}`Equation 3.19 <intext_eq_319_md>` for more details.

:::

(eq-322)=
:::{admonition} Minimized kriging estimate variance

$$
    \min \sigma_E^2(\mathbf{u}) = \min \mathbb{V}[\mathbf{z}^*(\mathbf{u}) - Z(\mathbf{u})]
$$

See {ref}`Equation 3.22 <intext_eq_322_md>` for more details.

:::

(eq-323)=
:::{admonition} Kriging linear matrix

\begin{equation*}
    \mathbf{\Gamma}_{n \times n} = 
        \begin{bmatrix}
            \gamma(\mathbf{u}_1 - \mathbf{u}_1) & \gamma(\mathbf{u}_1 - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_1 - \mathbf{u}_n) \\
            \gamma(\mathbf{u}_2 - \mathbf{u}_1) & \gamma(\mathbf{u}_2 - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_2 - \mathbf{u}_n) \\
            \vdots & \vdots & \ddots & \vdots \\
            \gamma(\mathbf{u}_n - \mathbf{u}_1) & \gamma(\mathbf{u}_n - \mathbf{u}_2) & \cdots & \gamma(\mathbf{u}_n - \mathbf{u}_n)
        \end{bmatrix}
\end{equation*}

See {ref}`Equation 3.23 <intext_eq_323_md>` for more details.

:::

(eq-324a)=
:::{admonition} Kriging estimation variance

$$
    \sigma_E^2(\mathbf{u}) = 
        - \sum_{i=1}^n \sum_{j=1}^n 
    \lambda_i \lambda_j 
    \gamma(\mathbf{u}_i - \mathbf{u}_j) +
        2 \sum_{j=1}^n \lambda_i \gamma(\mathbf{u}_i - \mathbf{u})
$$

See {ref}`Equation 3.24a <intext_eq_324a_md>` for more details.

:::

(eq-326b)=
:::{admonition} Simplified Lagrangian multiplier

$$
    \sum_{j=1}^n \lambda_j \gamma(\mathbf{u}_j - \mathbf{u}_b) + \mu = \gamma(\mathbf{u}_b - \mathbf{u})
$$

See {ref}`Equation 3.26b <intext_eq_326b_md>` for more details.

:::

(eq-326c)=
:::{admonition} Kriging optimization constraint

$$
    \frac{\partial \mathcal{L}}{\partial \mu} = 2 \left(1 - \sum_{j=1}^n \lambda_j\right) = 0
$$

See {ref}`Equation 3.26c <intext_eq_326c_md>` for more details.

:::

(eq-327)=
:::{admonition} Kriging system of equations required for solving for $\lambda$

\begin{equation*}
    \begin{cases}
        \sum\limits_{j=1}^n \lambda_j(\mathbf{u})\gamma(\mathbf{u}_j - \mathbf{u}_i) + \mu = \gamma(\mathbf{u}_i - \mathbf{u}) & \text{for } i = 1, \ldots, n \\
        \sum\limits_{j=1}^n \lambda_j(\mathbf{u}) = 1
    \end{cases}
    \tag{3.27}
\end{equation*}

See {ref}`Equation 3.27 <intext_eq_327_md>` for more details.

:::