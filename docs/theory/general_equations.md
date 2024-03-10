# Other summary statistics

![ text ](../images/example_indexing.jpg)


Length binning/quantization:

Given a distribution of length measurements:

$$
L = [ L_{1}, L_{2} , \dots , L_{n}]
$$

Given a distribution of age measurements:

$$
\delta = [ \delta_{1}, \delta_{2} , \dots , \delta_{n} ]
$$

where $n$ is $n$<sup>th</sup> measurement. 

Given the case where both $L_{\mathrm{cm}}$ and $\delta$ are binned:

$$
\varphi = \{ \left[ \iota_{a,1}, \iota_{b,1} \right] , \left[ \iota_{a,2}, \iota_{b,2} \right] , \dots , \left[ \iota_{a,m}, \iota_{b,m} \right] \}
$$

where $\varphi$ represents the binned intervals of a given variable, $\iota_{a,m}$ and $\iota_{b,m}$ are the infimum and supremum of the $m$<sup>th</sup> interval, respectively.
Both age and length measurements are then quantized:

$$
Q(x)_{n} = 
    \begin{cases} 
    \varphi_{m} & \text{if } \iota_{a,m} \leq x_{n} \leq \iota_{b,m} \text{ for some } m \text{ where } 1 \leq m \leq m_{\mathrm{max}} \\
    \text{undefined} & \text{otherwise}
\end{cases}
$$

where $Q(x)_{n}$ represents the quantization function applied to the $n$<sup>th</sup> values of either $L$ or $\delta$, $k$ denotes the indexed $\varphi$ interval.
The quantized age ($\alpha$) and length ($\ell$) measurements are then used to index biomass calculations.

Length-weight regression

$$L_{n} = \left( \frac{w_{n}}{b} \right) ^{\frac{1}{a}}$$