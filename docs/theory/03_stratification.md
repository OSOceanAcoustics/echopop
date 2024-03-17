# Stratification and apportioning of biological estimates


## Stratification
In practice, the acoustic quantities and biological estimates discussed in previous sections can vary depending on geospatial variation of the biological aggregations themselves. For example, the size and age of fish can vary depending on the survey location, as well as the sex of the fish. Therefore, the acoustic measurements and biological samples are typically stratified to account for these variations, and the biomass density is a function of the stratum (which depends on the geospatial locations) and sex, i.e.

$$
\rho_{w; s} = \rho^i_{w; s}(x,y) = \rho_A(x,y) (\mathbf{L}^i_s)^\top \mathbf{w}^i_s,
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
\rho_{w; s,\alpha} = \rho^i_{w; s,\alpha}(x,y) = \rho_A(x,y) (\mathbf{L}^i_{s,\alpha})^\top \mathbf{w}^i_{s,\alpha},
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



## Apportioning of kriged biomass

In Echopop, the kringing procedure interpolates the biomass density $\rho_w$ derived based on $NASC$ to finer grids where acoustic data are not collected. Let $\hat{\rho}_w$ be the kriged biomass density. The biomass of fish of sex $s$, length $\ell$, and age $\alpha$ at the kriging grid $(x_k, y_k)$ can be obtained by

$$
W_{s,\alpha,\ell}(x_k, y_k) = A(x_k, y_k) \; L^i_{s, \alpha, \ell} \; \hat{\rho}_w(x_k, y_k),
$$

where $A(x_k, y_k)$ and $\hat{\rho}_w(x_k, y_k)$ are the area and biomass density at the kriging grid $(x_k, y_k)$.

Note that the kriging grids can only be stratified with the INPFC straficiation (see below) based on the grid location, which also determins the stratum $i$ of the grid. The grid stratum in turn determins the number frequency of fish length $L^i_{s, \alpha, \ell}$ used in the apportioning.




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





## Hake survey specifics

### Stratification
For Pacific hake, two types of stratifications are used:

- **INPFC**: Stratification set by the International North Pacific Fisheries Commission (INFPC) that is based solely on latitude. The US-Canada bienniel hake survey region encompasses 6 strata.
- **KS**: Stratification determined based on the Kolomogorov-Smirnov test for differences of the fish length distributions across survey hauls.

### Haul sample measurements

After each haul, the hake samples are processed at two stations:

- Station 1: the length, weight, sex, and age of each fish are measured, and additional tissue samples are collected
- Station 2: only the length, weight, and sex of each fish are measured
