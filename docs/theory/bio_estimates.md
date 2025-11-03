(bio-estimates)=
# Biological estimates

## Biomass estimates

We can obtain an estimate of biomass density ($\rho_\text{B}$, kg nmi<sup>-2</sup>) by multiplying the [areal number density](./acoustics.md#eq-nasc-def) ($\rho_\text{A}$, animals nmi<sup>-2</sup>) of animals by the average weight ($\left< w \right>$, units: kg)

$$
\rho_\text{B} = \rho_\text{A} \left< w \right>,
\tag{1}
$$

where $\left< w \right>$ is the average weight: 

$$
\left< w \right> = \frac{\sum\limits_{j=1}^N w_j}{N}.
\tag{2}
\label{eq:average_weight}
$$

Here $w_j$ is the weight of fish $j$, and $N$ is the total number of fish samples. Many fisheries surveys discretize the observed continuous length ($L$) distributions, which allows $\eqref{eq:average_weight}$ to be redefined as: 

$$
\left< w \right> = \tilde{\mathbf{L}}^{\mathsf{T}} \mathbf{w},
\tag{3}
$$

where $\textbf{w}$ is the vector of summed $w$ in length bin $\ell$ and $\tilde{\mathbf{L}}$ the vector representing the normalized number frequency $\tilde{L}_\ell$ of fish samples in $\ell$:

\begin{align*}
    \tilde{\mathbf{L}} =
    \left[
        \begin{split}
            \tilde{L}&_1 \\ \tilde{L}&_2 \\ \tilde{L}&_3 \\ \vdots&
        \end{split}
    \right],
    \tag{4a}
\end{align*}

where if $L_\ell$ is nominal frequency of fish counts per $\ell$: 

\begin{align*}
    \begin{split}
        \tilde{L}_\ell = 
        \frac{L_\ell}{\sum\limits_{\ell} L_\ell},
    \end{split}
    \tag{4b}    
    \\[2ex]
    \begin{split}
        \sum\limits_\ell \tilde{L}_\ell = 1. 
    \end{split}
    \tag{4c}
\end{align*}

The weight vector, $\textbf{w}$, similarly represents the summed $w$ for each $\ell$:

\begin{align*}
    \mathbf{w} =
    \left[
        \begin{split}
            w&_1 \\ w&_2 \\ w&_3 \\ \vdots&
        \end{split}
    \right]
    \tag{5}
\end{align*}

The $w_\ell$ values are estimated by either summing the weights of fish belonging to each $\ell$, or fitting a log-linear regression to specimen length and weight measurements derived from trawl samples:

<a id="eq:ts_l_regression"></a>

$$
\log_{10}(\hat{w}) =
\log_{10}(\hat{a}) + \hat{b} \log_{10}(L),
\tag{6}
$$

where $\hat{w}$, $\hat{a}$, and $\hat{b}$ are the weight, $y$-intercept, and slope estimates fitted using ordinary least squares (OLS). If $L_{\ell}^*$ is the representive length for bin $\ell$, then $w_\ell$ can be calculated by:

$$
w_\ell =
\left[
    10^{\hat{a}} {(L_{\ell}^{*})}^{\hat{b}}
\right]
L_\ell.
\tag{7}
$$

With the above quantities, the biomass ($B$, kg) can then be estimated by:

$$
B = \rho_B A = \rho_A \left< w \right> A,
\tag{8}
$$

where $A$ is the unit area (nmi<sup>2</sup>) associated with the density measure.