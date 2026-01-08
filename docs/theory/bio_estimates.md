(bio-estimates)=
# Biological estimates

## Biomass estimates

We can obtain an estimate of biomass density ($\rho_\text{B}$, kg nmi<sup>-2</sup>) by multiplying the areal number density ({ref}`Eq. 1.15 <eq-115>`) of animals by the average weight $\left< w \right>$ (kg):

<a id="intext_eq_21"></a>

(intext_eq_21_md)=
\begin{align*}
    \rho_\text{B} = \rho_\text{A} \left< w \right>,
    \tag{2.1a}
    \\[2ex]
    \left< w \right> = \frac{\sum\limits_{j=1}^N w_j}{N}.
    \tag{2.1b}
\end{align*}

Here $w_j$ is the weight of fish $j$ and $N$ is the total number of fish samples. Many fisheries surveys discretize the observed continuous length ($L$) distributions, which allows {ref}`Eq. (2.1b) <eq-21b>` to instead be expressed as:

<a id="intext_eq_22"></a>

(intext_eq_22_md)=
$$
    \left< w \right> = \tilde{\mathbf{L}}^{\mathsf{T}} \mathbf{w},
    \tag{2.2}
$$

where $\textbf{w}$ is the vector of summed $w$ in length bin $\ell$ and $\tilde{\mathbf{L}}$ the vector representing the normalized number frequency $\tilde{L}_\ell$ of fish samples in $\ell$:

<a id="intext_eq_23"></a>

(intext_eq_23_md)=
\begin{align*}
    \tilde{\mathbf{L}} =
    \left[
        \begin{split}
            \tilde{L}&_1 \\ \tilde{L}&_2 \\ \tilde{L}&_3 \\ \vdots&
        \end{split}
    \right].
    \tag{2.3}
\end{align*}

Let $L_\ell$ denote the nominal frequency (or count) of fish per $\ell$. The normalized frequency is then:

\begin{align*}
    \begin{split}
        \tilde{L}_\ell = 
        \frac{L_\ell}{\sum\limits_{\ell} L_\ell},
    \end{split}
    \tag{2.4a}    
    \\[2ex]
    \begin{split}
        \sum\limits_\ell \tilde{L}_\ell = 1. 
    \end{split}
    \tag{2.4b}
\end{align*}

The weight vector, $\textbf{w}$, similarly represents the summed $w$ for each $\ell$:

<a id="intext_eq_25"></a>

(intext_eq_25_md)=
\begin{align*}
    \mathbf{w} =
    \left[
        \begin{split}
            w&_1 \\ w&_2 \\ w&_3 \\ \vdots&
        \end{split}
    \right].
    \tag{2.5}
\end{align*}

Values of $w_\ell$ are estimated by either summing the weights of fish belonging to each $\ell$, or fitting a log-linear regression to specimen length and weight measurements derived from trawl samples:

<a id="intext_eq_26"></a>

(intext_eq_26_md)=
$$
    \log_{10}(\hat{w}) =
        \log_{10}(\hat{a}) + \hat{b} \log_{10}(L),
    \tag{2.6}
$$

where $\hat{w}$, $\hat{a}$, and $\hat{b}$ are the weight, $y$-intercept, and slope estimates fitted using ordinary least squares (OLS). Let $L_{\ell}^*$ denote the representive length for bin $\ell$. Then the weight for bin $\ell$ is:

<a id="intext_eq_27"></a>

(intext_eq_27_md)=
$$
    \mathcal{W}(\ell) =
        \left[
            10^{\hat{a}} {(L_{\ell}^{*})}^{\hat{b}}
        \right]
        L_\ell.
    \tag{2.7}
$$

With the above quantities, the biomass ($B$, kg) can then be estimated by:

$$
    B = \rho_\text{B} A = \rho_\text{A} \left< w \right> A,
    \tag{2.8}
$$

where $A$ is the unit area (nmi<sup>2</sup>) associated with the areal density estimate ({ref}`Eq. 1.15 <eq-115>`).