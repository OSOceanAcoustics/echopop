(semivariogram_eq)=
# Semivariogram equations


Spherical:

$$
\gamma(h) = C(0) \left[1.5\frac{h}{L} - 0.5\left(\frac{h}{L}\right)^3\right] + \gamma(0)
$$

Exponential:

$$
\gamma(h) = C(0) (1-e^{-h/L}) + \gamma(0)
$$

Gaussian:

$$
\gamma(h) = C(0) (1-e^{-(h/L)^2}) + \gamma(0)
$$

Linear:

$$
\gamma(h) = C(0) h + \gamma(0)
$$

Sinc:

$$
\gamma(h) = C(0) \left[1-\sin(bh)\right] + \gamma(0)
$$

Composite exponential-cosine:

$$
\gamma(h) = C(0) \left[1 - \cos(b h) e^{-h/L} \right] + \gamma(0)
$$

Composite exponential-cosine with enhanced semivariance

$$
\gamma(h) = C(0) \left[1 +  \cos(b h)e^{-h/L}\right] + \gamma(0)
$$

Composite Gaussian-cosine 

$$
\gamma(h) = C(0) \left[1- \cos(b h)e^{-(h/L)^2}\right] + \gamma(0)
$$

Bessel

$$
\gamma(h) = C(0) \left( 1 - J_0(bh)\right) + \gamma(0)
$$

Composite Bessel-exponential:

$$
\gamma(h) = C(0) \left[1- J_0(b h)e^{-h/L}\right] + \gamma(0)
$$

Composite Gaussian-Bessel:

$$
\gamma(h) = C(0) \left[1- J_0(b h)e^{-(h/L)^2}\right] + \gamma(0)
$$

Composite Gaussian-linear:

$$
\gamma(h) = C(0) \left[1- (1- b h) e^{-(h/L)^2}\right] + \gamma(0)
$$

Generalized Bessel-exponential composite model:

$$
\gamma(h) = C(0) \left[1 - J_0(bh)e^{-(h/L)^p}\right] + \gamma(0)
$$
