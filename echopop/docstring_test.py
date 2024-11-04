        """
        Plotting method for visualizing results

        Note
        ----------
        *See `Plotting parameterization` section for more details*

        Parameters
        ----------
        kind: Literal["age_length_distribution", "mesh", "transect"]
            The 'kind' of plot and dataset that will be visualized. Possible plots include:
            
            - `age_length_distribution`: two-dimensional heatmap of population estimates
                as a function of age (x-axis) and length (y-axis).    
                                        
            - `mesh`:  map displaying kriged estimates distributed across the spatial kriging
                mesh.            
                
            - `transect`: bubbleplot of various along-transect population estimates.     
                           
        variable: str
            The variable used for plotting. *See the `Plotting Variables` section for more details*.
        plot_parameters: Dict[str, Any]
            A dictionary comprising various plotting parameters. *See the `Plotting parameters`
            for more details*.
        plot_type: Optional[Literal["heatmap", "hexbin", "scatter", "pcolormesh"]]
            The type of plot. Options are specific to each `kind`:
            
            - `age_length_distribution`: ["heatmap"]                
            - `mesh`: ["hexbin", "scatter", "pcolormesh"]                
            - `transect`: ["scatter"]

        Other Parameters
        ----------------
        Valid variables depend on the input for `kind`:

        "age_length_distribution":
        
            - **"abundance"**: Animal abundance.
            - **"biomass"**: Animal biomass.
            
        "mesh":
        
            - **"biomass"**: Kriged animal biomass.
            - **"kriged_mean"**: Kriged animal biomass density.
            - **"kriged_variance"**: Kriging variance.
            - **"sample_cv"**: Mesh node coefficient of variation.
            - **"sample_variance"**: Sample variance.
            
        "transect":
        
            - **"abundance"**: Animal abundance.
            - **"abundance_female"**/**"abundance_male"**: Sexed animal abundance.
            - **"biomass"**: Animal biomass.
            - **"biomass_female"**/**"biomass_male"**: Sexed animal biomass.
            - **"biomass_density"**: Animal biomass density.
            - **"biomass_density_female"**/**"biomass_density_male"**: Sexed animal biomass density.
            - **"nasc"**: Nautical area scattering coefficient (NASC).
            - **"number_density"**: Animal number density.
            - **"number_density_female"**/**"number_density_male"**: Sexed animal number density.
        
        The format for `plotting_parameters` is:

        axis_limits: Optional[Dict[str, Any]]
        
            A dictionary that contains limits for the x- and y-axes. This should contain two
            nested dictionaries for `x` and `y`:
            
            - `xmin`/`ymin` (float): Minimum x- and y-axis values.            
            - `xmax`/`ymax` (float): Maximum x- and y-axis values.            
            - `left`/`right` (float): The left- and right-most axis values. This should not be
            included if `xmin`/`xmax`/`ymin`/`ymax` are included.
            
        geo_config: Optional[Dict[str, Any]]
            A dictionary that contains three keyword arguments:
            
            - `coastline` (`cartopy.feature.Feature`): A `cartopy.feature.Feature` object that
            includes land and/or coastline information used for geospatial plots.            
            - `init` (str): The initial geospatial projection code (e.g. "EPSG:4326").            
            - `plot_projection` (`cartopy.crs.Projection`): A `cartopy.crs.Projection`
            projection object used for plotting geospatial data.
            
        grid_heatmap: Optional[bool]
            An optional boolean value for adding correctly spaced gridlines to the biological
            heatmap. It is otherwise not used for the spatial plots.
        sex: Optional[Literal["all", "female", "male"]]
            An optional argument for defining which fish sex to plot for the biological
            heatmap. It is otherwise not used for the spatial plots.
        log_base: Optional[float]
            Base for rescaling the plot colormap via logarithmic transformation.
        cmap: Optional[str]
            Plotting colormap (e.g. "viridis").
        vmin, vmax: Optional[float]
            The data range used to rescale the colormap in linear or logarithmic space.
        **kwargs: Dict[str, Any]
            Plotting functions can accept additional arguments used within
            :func:`matplotlib.pyplot.plot`
        
        See Also
        ----------------
        :func:`matplotlib.pyplot.plot`
            For more details on `matplotlib` plotting function keyword arguments. Keep in mind that
            these should be used in caution as they may interfere with the primary variables
            defined for `plotting_parameters`.

        """