(stratification)=
# Stratification to obtain biomass density


(stratification-intro)=
## Stratification
In practice, the acoustic quantities and biological estimates discussed in previous sections can vary depending on geospatial variation of the biological aggregations themselves. For example, the size and age of fish can vary depending on the survey location, as well as the sex of the fish. Therefore, the acoustic measurements and biological samples are typically stratified to account for these variations, and the biomass density is a function of the stratum (which depends on the geospatial locations) and sex, i.e.

$$
\rho_{B; s} = \rho^i_{B; s}(x,y) = \rho_A(x,y) (\mathbf{L}^i_s)^\top \mathbf{w}^i_s,
$$

where $i$ is the stratum, $\rho_A(x,y)$ is the nautical areal number density at location $(x, y)$, and $\mathbf{L}^i_s$ and $\mathbf{w}^i_s$ are the vectors characterizing the number frequency of fish length and the corresponding weight in stratum $i$, for fish of sex $s$:

$$
\mathbf{L}^i_s = \begin{bmatrix}
L^i_{s,1} \\
L^i_{s,2} \\
L^i_{s,3} \\
\vdots
\end{bmatrix},
$$

and

$$
\mathbf{w}^i_s = \begin{bmatrix}
w^i_{s,1} \\
w^i_{s,2} \\
w^i_{s,3} \\
\vdots
\end{bmatrix}.
$$


Note that the number frequency of fish here is normalized across all length bins and sex, i.e., 

$$
\sum_{s,\ell} L^i_{s,\ell} = 1
$$


### Including age data
In the case when fish age is measured and binned, the biomass density is a function of the stratum (which depends on the geospatial locations), sex, and age:

$$
\rho_{B; s,\alpha} = \rho^i_{B; s,\alpha}(x,y) = \rho_A(x,y) (\mathbf{L}^i_{s,\alpha})^\top \mathbf{w}^i_{s,\alpha},
$$

where $\alpha$ is the age bin,

$$
\mathbf{L}^i_{s,\alpha} = \begin{bmatrix}
L^i_{s,\alpha,1} \\
L^i_{s,\alpha,2} \\
L^i_{s,\alpha,3} \\
\vdots
\end{bmatrix},
$$

and 

$$
\mathbf{w}^i_{s,\alpha} = \begin{bmatrix}
w^i_{s,\alpha,1} \\
w^i_{s,\alpha,2} \\
w^i_{s,\alpha,3} \\
\vdots
\end{bmatrix}.
$$

All of $L^i_{s,\alpha,\ell}$ and $w^i_{s,\alpha,\ell}$ vary depending on the stratum $i$, the fish sex $s$, and the age bin $\alpha$.


Note that the number frequency of fish length here is normalized across all age bins, length bins, and sex within a stratum, i.e.

$$
\sum_{s,\ell,\alpha} L^i_{s,\alpha,\ell} = 1
$$



<!-- ## Jolly-Hampton (1990) stratified sampling 

Mean density for stratum $i$:

$$
\hat{ \rho }_{A,B}^{ i } = 
    \frac{1}{ n^{ i } }
    \sum\limits_{i=0}^{n^{i} } w^{i,j} \hat{ \rho }_{A,B}^{ i,j,k}
$$

where $w^{i,j}$ is the transect weight calculated via:

$$
w^{i,j} = \frac{
    d(x,y)^{i,j}
}{
    \frac{1}{n^{i}}
    \sum\limits_{j=1}^{n^{i}} d(x,y)^{i,j}
}
$$

where $d(x,y)^{i,j}$ is the transect length of $n^{i}$ transects within each stratum. 
This procedure is then repeated by using the different areas ($A^{i}$) of each stratum to
weight the final $\hat{ \rho }_{A,B}$ estimate:

$$
\hat{ \rho }_{A,B} =
    \frac{
        \sum\limits_{i} A_{i} \hat{ \rho_{A,B}^{ i } }
    }{
        \sum\limits_{i} A_{i}
    }
$$ -->


## Stratification schemes used in the hake survey
For Pacific hake, two types of stratifications are used:

- **INPFC**: Stratification set by the International North Pacific Fisheries Commission (INFPC) that is based solely on latitude. The US-Canada bienniel hake survey region encompasses 6 strata.
- **KS**: Stratification determined based on the Kolomogorov-Smirnov test for differences of the fish length distributions across survey hauls.
