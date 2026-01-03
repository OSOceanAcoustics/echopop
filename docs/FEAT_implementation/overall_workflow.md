# Overall data workflow

The Echopop data workflow consists of 5 primary components:
- Data ingestion
    - Set up paths
    - Ingest NASC, region, haul data
    - Ingest biological data
    - Load stratification schemes
    - Load kriging parameters
- Transect interval-based data processing
    - Stratify ingested data
    - Compute length-weight relationship
    - Compute number proportion
    - Compute weight proportion
- Biological estimates (transect interval-based)
    - Compute various biological quantities with following sequence: NASC ($\rho_A$) → abundance ($\rho_N$) → biomass density $\rho_B$ → biomass ($B$) (see figure below)
    - Distribute biological estimates across (sex, length, age) for each transect interval
- Kriging
    - Coordinate transformation 
    - Variogram analysis and fitting
    - Perform kriging on biomass density
    - Derive/back-calculate other biological estimates: kriged biomass density ($\hat{\rho}_B$) → kriged biomass ($\hat{B}$) → kriged abundance ($\hat{\rho}_N$) → kriged NASC ($\hat{\rho}_A$) (see figure below)
    - Distribute kriged biological estimates across (sex, length, age)
- Stratified analysis
    - Perform Jolly-Hampton analysis using transect-based data and kriged data

The documentation is organized with roughly the same order to explain each component in detail.

![](../assets/echopop_data_flow.png)