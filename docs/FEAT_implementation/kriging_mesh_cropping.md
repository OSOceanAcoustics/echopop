(kriging-mesh-cropping)=

# Kriging mesh cropping

When cropping the kriging mesh, this is done using the `Kriging` class where:

```python
from echopop.geostatistics import Kriging
import echopop.workflows.nwfsc_feat as feat

# Define the requisite kriging parameters
KRIGING_PARAMETERS = {
    "search_radius": best_fit_parameters["correlation_range"] * 3,
    "aspect_ratio": 0.001,
    "k_min": 3,
    "k_max": 10,
}

# Define the requisite variogram parameters and arguments
VARIOGRAM_PARAMETERS = {"model": ["exponential", "bessel"], **best_fit_parameters}

# Create Kriging instance
kriger = Kriging(
    mesh=df_mesh,
    kriging_params=KRIGING_PARAMETERS,
    variogram_params=VARIOGRAM_PARAMETERS,
    coordinate_names=("x", "y"),
)
```

Here, `df_mesh` corresponds to a `pandas.DataFrame` with a series of coordinates that are used 
for kriging interpolation. The `Kriging` class defaults to `extrapolation = True` when using the `Kriging.krige` method:

```python
# Execute ordinary kriging
df_kriged_results = kriger.krige(
    transects=df_nasc,
    variable="biomass_density",
    extrapolate=True,
    default_mesh_cell_area=6.25,
)
```

where `df_nasc` corresponds to a `pandas.DataFrame` comprising the transect data. However, when `extrapolate=False`, the default method for cropping the mesh uses a convex hull approach. Prior to running `Kriging.krige`, the `Kriging.crop_mesh` method can be applied:

```python
# Crop mesh using the default approach
kriger.crop_mesh(
    transects=df_nasc,
    num_nearest_transects = 3,
    mesh_buffer_distance = 2.5,
)
```

## Extending cropping methods

Custom cropping methods can be used instead of the default convex hull approach by using the `Kriging.crop_mesh(..., crop_function)` argument. This accepts any `Callable` function. For instance, the cropping method from EchoPro is implemented via the `transect_ends_crop()` function:

```python
# Crop mesh using the transect-ends approach
kriger.crop_mesh(
    transects=df_nasc,
    latitude_resolution=1.25 / 60.0,
    mesh_buffer_distance = 2.5,
    transect_mesh_region_function=feat.parameters.transect_mesh_region_2019,
)
```

Here, the <b>only</b> required argument in `Kriging.crop_mesh` is `crop_function`. The keyword arguments for the defined function can then be supplied to `Kriging.crop_mesh`, such as these keyword arguments required by `transect_ends_crop`:

- `transects`: `pandas.DataFrame` comprising the transect data with the same coordinates (projection/transformation) as the mesh grid
- `mesh`: `pandas.DataFrame` comprising the mesh grid coordinates. Note that this was not specified in the code snippet above. This is because the mesh grid is already stored internally within the `Kriging`-class object.
- `latitude_resolution`: The latitudinal resolution (in degrees) used for the interpolation. This determines the spacing between interpolation points and affects the precision of the boundary detection.
- `transect_mesh_region_function`: This is a `Callable` function that sorts and maps each transect number to a discretized region in the mesh grid. 

## Transect mesh region functions

These functions mirror the `transect_region_def_*.m` files used by EchoPro. These are tailored specifically to each survey, but can otherwise all be imported directly from `echopop.workflows.nwfsc_feat.parameters`. As an example, the `transect_mesh_regiuon_function` for the 2019 survey is:

```python
from typing import List, Tuple

def transect_mesh_region_2019(
    region: np.number,
) -> Tuple[np.number, np.number, List[np.number], List[np.number]]:

    # Initialize variables
    transect_start = None
    transect_end = None
    transect_lower_bound = []  # W/S
    transect_upper_bound = []  # E/N

    # Region 1: parallel transects to latitudes from south of SCB to west of Haida Gwaii
    if region == 1:
        # ---- Southern-most transect
        transect_start = 1
        # ---- Northern-most transect
        transect_end = 119
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)]
    # Region 2: transects parallel to longitudes north of Haida Gwaii
    elif region == 2:
        # ---- Western-most transect
        transect_start = 121
        # ---- Eastern-most transect
        transect_end = 127
        # ---- Southern boundary
        transect_lower_bound = [i + 0.6 for i in range(transect_start, transect_end + 1)]
        # ---- Northern boundary
        transect_upper_bound = [i + 0.9 for i in range(transect_start, transect_end + 1)]
    # Region 3: parallel transects to latitudes west of Haida Gwaii
    else:
        # ---- Southern-most transect
        transect_start = 129
        # ---- Northern-most transect
        transect_end = 145
        # ---- Western boundary
        transect_lower_bound = [i + 0.1 for i in range(transect_start, transect_end + 1)]
        # ---- Eastern boundary
        transect_upper_bound = [i + 0.4 for i in range(transect_start, transect_end + 1)]

    return transect_start, transect_end, transect_lower_bound, transect_upper_bound
```

This returns a `tuple` that includes the first and last transects of each region, and those that make up the lower- and uppermost boundaries (i.e., either east-west or north-south, depending on the specification).

## Other custom methods for kriging

Along with being able to use custom cropping functions, user-created adaptive search strategies can be registered for use by the kriging algorithm. These adaptive methods are how the algorithm searches for the nearest neighbors and resolve weights assigned to values at each coordinate. The `western_boundary_search_strategy` function represents the approached used by EchoPro and can be added to the `Kriging`-class object via:

```python
kriger.register_search_strategy("FEAT_strategy", feat.western_boundary_search_strategy)
```

Once registered, it can be applied via:

```python
# Define the required keyword arguments for 'FEAT_strategy'
# ---- Only `transect_western_extents` is needed for this particular function since the
# `kriging_mesh` and `coordinate_names` arguments are inherited from the class instance
FEAT_STRATEGY_KWARGS = {
    "western_extent": transect_western_extents,
}

# Apply ordinary kriging
df_kriged_results = kriger.krige(
    transects=df_nasc,
    variable="biomass_density",
    extrapolate=False,
    default_mesh_cell_area=6.25,
    adaptive_search_strategy="FEAT_strategy",
    custom_search_kwargs=FEAT_STRATEGY_KWARGS,
)
```

Similar to the `Kriging.crop_mesh` method, keyword arguments specific to the custom search strategy can be provided; however, they are instead contained within the `custom_search_kwargs` dictionary. Otherwise, custom functions can utilize a variety of internal variables (e.g. `k_min`, `oos_weights`) depending on user needs. In this case, the only external keyword argument required for `"FEAT_strategy"` is `transect_western_extents`. This is a `pandas.DataFrame` generated by the `get_survey_western_extents` function:

```python
transect_western_extents = feat.get_survey_western_extents(
    transects=df_nasc, coordinate_names=("x", "y"), latitude_threshold=51.0
)
```

This gets the westernmost extent of each transect that is used to constrain the adaptive nearest neighbors search strategy incorporated into the kriging interpolation algorithm. 