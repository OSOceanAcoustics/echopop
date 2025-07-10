(kriging-eq)=
# Kriging equations

Intrinsic model 

$$
\mathbb{E} \left[ ( z( \boldsymbol{ \mathrm{ x } } ) - m )( z( \boldsymbol{ \mathrm{ x } } + \boldsymbol{ \mathrm{ h } } ) - m ) \right] = C(h)
$$

$$
h = \lVert \boldsymbol{ \mathrm{ x } } - ( \boldsymbol{ \mathrm{ x } } + \boldsymbol{ \mathrm{ h } } ) \rVert =
\sqrt{
    ( x_{1} - x'_{1} ) ^ 2 + 
    ( x_{2} - x'_{2} ) ^ 2 +
    ( x_{3} - x'_{3} ) ^ 2
}
$$

Expression of the Intrisic model via semi-variogram:

$$
\gamma(h) = \frac{1}{2} \mathbb{E}\left[ z(\boldsymbol{ \mathrm{ x } } ) - z( \boldsymbol{ \mathrm{ x } } + \boldsymbol{ \mathrm{ h } } ) ^ 2 \right]
= C(0) - C(h)
$$

Weighted averaging: 

$$
\hat{ z }( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) =
    \sum\limits_{\alpha}^{N} \lambda_{\alpha} z(\boldsymbol{ \mathrm{ x } }_\alpha)
$$

where $\lambda_{\alpha}$ is the weighting coefficient. This quantity can be estimated from the variance ($\mathbb{V}$):

$$
\begin{array}{l}
    \mathbb{V}( z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) ) =  
        \mathbb{E} \{ \left[  z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) - \hat{ z }( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) \right ] ^ 2 \} \\
    \phantom{\mathrm{ Var( z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) ) }} = 
        C(0) - 2 \sum\limits_{\alpha} \lambda_{\alpha} 
        C( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \kappa } \rVert )
        + \sum\limits_{ \alpha } \sum\limits_{ \beta } \lambda_{ \alpha } \lambda_{ \beta }
        C( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \beta } \rVert )
\end{array}
$$

$$
\sum\limits_{\alpha}^{N} \lambda_{\alpha} = 1
$$

Differentiation with respect to $\lambda_\alpha$ to find the predicted value while also minimizing the variance:

$$
\sum\limits_{\beta}^{ N } \lambda_{ \beta } 
    C_{n}( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \beta } \rVert )
    - \mu = 
    C_{n}( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \kappa } \rVert )
$$

$$
\sum\limits_{\beta}^{N} \lambda_{\beta} = 1
$$

This can also be expressed as a system of linear equations that can be solved to determine $\lambda_{\alpha}$ and 
$\lambda_{\beta}:$

$$
\begin{bmatrix}
    \gamma(d(u_1,u_1)) & \gamma(d(u_1,u_2)) & \cdots & \gamma(d(u_1,u_n)) & 1 \\
    \gamma(d(u_2,u_1)) & \gamma(d(u_2,u_2)) & \cdots & \gamma(d(u_2,u_n)) & 1 \\
    \vdots & \vdots & \ddots & \vdots & \vdots \\
    \gamma(d(u_n,u_1)) & \gamma(d(u_n,u_2)) & \cdots & \gamma(d(u_n,u_n)) & 1 \\
    1 & 1 & \cdots & 1 & 0
    \end{bmatrix}
    \begin{bmatrix}
    \lambda_1 \\
    \lambda_2 \\
    \vdots \\
    \lambda_n \\
    \mu
    \end{bmatrix}
    =
    \begin{bmatrix}
    \gamma(d(u,u_1)) \\
    \gamma(d(u,u_2)) \\
    \vdots \\
    \gamma(d(u,u_n)) \\
    1
\end{bmatrix}
$$

where $\mu$ is the Lagrangian coefficient:

$$
\mu = \frac{
    1 - 1^{T} \lambda 
}{
    \mathbf{1}^{T} K^{-1} \boldsymbol{ \mathrm{ k } }
}
$$

where $K$ is the $n~\mathrm{x}~n$ covariance matrix, $\boldsymbol{ \mathrm{ k } }$ is the 
$n~\mathrm{x}~1$ vector of covariance estimates between sampled and unsampled locations, 
$1$ corresponds to a vector of ones, and $\mathbf{1}^{T}$ is the transpose of the vector 
of ones. Thus, the above system of linear equations for $\lambda$ can be solved via:

$$
\lambda = K^{-1} (\mathbf{k} - \mu)
$$

This can also be expressed directly in terms of the semivariogram:

$$
\sum\limits_{\beta}^{ N } \lambda_{ \beta } 
    \gamma_{n}( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \beta } \rVert )
    + \mu = 
    \gamma_{n}( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } \boldsymbol{ \mathrm{ x } }_{ \kappa } \rVert )
$$

$$
\sum\limits_{\beta}^{N} \lambda_{\beta} = 1
$$

Once the $\lambda_{\beta}$ and $\mu$ estimates have been obtained, the kriged prediction variance can be 
estimated via:

$$
\begin{array}{l}
    \mathbb{V}( z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) ) =  
        \mathbb{E} \{ \left[  z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) - \hat{ z }( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) \right ] ^ 2 \} \\
    \phantom{\mathrm{ Var( z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) ) }} = 
        C(0) + \mu
        - \sum\limits_{\alpha} \lambda_{\alpha}
        C( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } - \boldsymbol{ \mathrm{ x } }_{ \kappa } \rVert ) \\
    \phantom{\mathrm{ Var( z( \boldsymbol{ \mathrm{ x } }_{ \kappa } ) ) }} = 
        \mu + \sum\limits_{\alpha} \lambda_{\alpha}
        \gamma( \lVert \boldsymbol{ \mathrm{ x } }_{ \alpha } - \boldsymbol{ \mathrm{ x } }_{ \kappa } \rVert )
        -\gamma(0)
\end{array}
$$

Ordinary kriging estimation equation:

$$
\hat{z}(u) = \sum_{i=1}^{n} \lambda_i z(u_i)
$$

$$
\sum\limits_{i=1}^{n} \lambda_{i} = 1
$$
