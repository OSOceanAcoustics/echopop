# Jolly and Hampton (1990) stratified sampling 

$i$ = stratum
$j$ = transect
$k$ = along-transect interval/segment
$s$ = sex
$\alpha$ = quantized age 
$\ell$ = quantized length

Mean density for the $i$<sup>th</sup> stratum.

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
$$