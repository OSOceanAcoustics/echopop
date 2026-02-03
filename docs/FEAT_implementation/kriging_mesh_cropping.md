(kriging-mesh-cropping)=

# Kriging mesh cropping and adaptive search strategies

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

where `df_nasc` corresponds to a `pandas.DataFrame` comprising the transect data. However, when `extrapolate=False`, the default method for cropping the mesh uses a convex hull approach. Prior to running `Kriging.krige`, the {meth}`Kriging.crop_mesh <echopop.geostatistics.Kriging.crop_mesh>` method can be applied:

(convex_hull_plot)=
```python
# Crop mesh using the default approach
kriger.crop_mesh(
    transects=df_nasc,
    num_nearest_transects = 3,
    mesh_buffer_distance = 2.5,
)
```

For reference, a convex hull is the smallest convex polygon that contains all the points in a given set. This boundary acts as a "hard limit" for the kriging interpolation, preventing the model from generating predictions in areas where no data exists. 

<iframe src="../_static/mesh_cropping_convex.html" width="650" height="400"></iframe>

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

<iframe src="../_static/mesh_cropping.html" width="650" height="400"></iframe>

```{tip} Comparing crop methods
How does this compare and/or contrast to the {ref}`convex hull approach <convex_hull_plot>`? There should be fairly reasonable agreement in the survey extents torward the centers of each transect mesh region.
```

## Transect mesh region functions

These functions mirror the `transect_region_def_*.m` files used by EchoPro. These are tailored specifically to each survey, but can otherwise all be imported directly from `echopop.workflows.nwfsc_feat.parameters`. For example, the 2019 survey `transect_mesh_region_function` is represented by the `transect_mesh_region_2019` function, which corresponds to the equivalent `transect_region_def_2019` function used by EchoPro (`MATLAB`):

::::{tab-set}
:::{tab-item} Echopop (Python)
:sync: tab1

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
:::
:::{tab-item} EchoPro (MATLAB)
:sync: tab2

```MATLAB
function [tx0,tx1,tx_out1,tx_out2]=transect_region_def_2019(region)

global para

tx0=[];
tx1=[];
tx_out1=[];
tx_out2=[];
if region == 1
%% region 1: paralell transects to latitudes from south of SCB to west of QC IS
    tx0=1;    % southern most transect number
    if para.proc.source == 1
        tx1=86;  % northern most transect number (US only)
    else
        tx1=119;  % northern most transect number (US & CAN)
    end
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% left (west) bound
    tx_l=[tx0:tx1]+0.1;
    %% right (east) bound
    tx_r=[tx0:tx1]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
elseif region == 2
%% region 2: transects paralell to longitudes north of QCI
    tx0=121;    % west most transect number
    tx1=127;    % east most transect number
%% specifies lower (south) and upper (north) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[tx0:tx1] + 0.6;
    tx_u=[tx0:tx1] + 0.9;
    tx_out1=tx_l;
    tx_out2=tx_u;
else
    %% region 3: paralell transects to latitudes west of QC IS
    tx0=129;    % northern most transect number
    tx1=145;    % southern most transect number
    %% specifies left (west) and right (east) region boundaries based on the transects
    %% #.1 = west end of transect
    %% #.4 = east end of transect
    %% #.6 = south end of transect
    %% #.9 = north end of transect
    tx_l=[tx0:tx1]+0.1;
    tx_r=[tx0:tx1]+0.4;
    tx_out1=tx_l;
    tx_out2=tx_r;
end

return
```
:::
::::

This returns a `tuple` that includes the first and last transects of each region, and those that make up the lower- and uppermost boundaries (i.e., either east-west or north-south, depending on the specification).

## Adaptive kriging search strategy

Along with being able to use custom cropping functions, user-created adaptive search strategies can be registered for use by the kriging algorithm. These adaptive methods are how the algorithm searches for the nearest neighbors and resolve weights assigned to values at each coordinate. The default search strategy is an {func}`uniform approach <echopop.geostatistics.uniform_search_strategy>` that {ref}`searches for the nearest neighbors <adaptive-search>` to constrain computational costs.

 The {func}`western_boundary_search_strategy <echopop.workflows.nwfsc_feat.functions.western_boundary_search_strategy>` function represents the approached used by EchoPro and can be added to the `Kriging`-class object via:

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

Similar to the `Kriging.crop_mesh` method, keyword arguments specific to the custom search strategy can be provided; however, they are instead contained within the `custom_search_kwargs` dictionary. Otherwise, custom functions can utilize a variety of internal variables (e.g. `k_min`, `oos_weights`) depending on user needs. In this case, the only external keyword argument required for `"FEAT_strategy"` is `transect_western_extents`. This is a `pandas.DataFrame` generated by the {func}`get_survey_western_extents <echopop.workflows.nwfsc_feat.functions.get_survey_western_extents>` function:

```python
transect_western_extents = feat.get_survey_western_extents(
    transects=df_nasc, coordinate_names=("x", "y"), latitude_threshold=51.0
)
```

This gets the westernmost extent of each transect that is used to constrain the adaptive nearest neighbors search strategy incorporated into the kriging interpolation algorithm. 

### Creating custom adaptive search strategies

Any adaptive nearest-neighbor search strategy can be registered to the `Kriging`-class object via `Kriging.register_search_strategy`. A custom function can be created tthat uses external arguments (e.g., `"western_extent"`) that can also reference intermediate variables internal to the `Kriging`-class object. These include:

- `sparse_radii`: An integer array comprising indices where there are fewer than ``k_min`` nearest neighbors.

- `valid_distances`: An integer array with the number of masked distance matrix values where extrapolation is required.

- `local_points`: An array with the sorted distances (from nearest to furthest) relative to each point.

- `distance_matrix_masked`: An array with the search-radius-masked nearest neighbor distances.

- `nearby_indices`: An integer array comprising indices of points that require extrapolation.

- `k_min`: The minimum number of nearest neighbors required for including values for kriging within
the search radius.

- `k_max`: The maximum number of nearest neighbors required for including values for kriging detected
within the search radius.

- `search_radius`: The adaptive search radius that identifies the *k*-nearest neighbors around each
georeferenced value that are subsequently kriged.

- `wr_indices`: Integer array containing the indices of within-radius (WR) (i.e. < ``k_max``) points.

- `oos_indices`: A template array based on the size of the data input and ``k_min`` that will contain indices
where extrapolation is required where there are fewer than ``k_min`` nearest neighbors.

- `oos_weights`: An array with weights that are applied to extrapolated values.

The only **required** arguments are `k_min` and `k_max`. The function **must** also a `tuple` comprising the following arrays:

- The *k*-nearest distances to each and every point.

- The within-radius indices that do not require any extrapolation.

- The  out-of-sample indices that will be extrapolated.

- The out-of-sample weights that are applied to all values that are extrapolated.

Comparisons between the default `uniform_search_strategy` and FEAT-specific `western_boundary_search_strategy` further illustrate potential places to edit the adaptive search strategy used by the ordinary kriging algorithm.

::::{tab-set}
:::{tab-item} Uniform
:sync: tab1

```python
def uniform_search_strategy(
    sparse_radii: np.ndarray[int],
    valid_distances: np.ndarray[int],
    local_points: np.ndarray[float],
    k_min: int,
    k_max: int,
    search_radius: float,
    wr_indices: np.ndarray[int],
    oos_indices: np.ndarray[np.number],
    oos_weights: np.ndarray[float],
    **kwargs,
) -> Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]:
    """
    Uniform extrapolation search strategy for finding (and weighting) k-th nearest points
    (relative to a reference coordinate) required for computing the lagged semivariogram in an
    adaptive approach

    Parameters
    ----------
    sparse_radii : |np.ndarray[int]|
        Indices where there are fewer than ``k_min`` nearest neighbors.
    valid_distances : |np.ndarray[int]|
        The number of masked distance matrix values where extrapolation is required.
    local_points : |np.ndarray[float]|
        An array with the sorted distances (from nearest to furthest) relative to each point.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.
    search_radius : float
        The adaptive search radius that identifies the :math:`k`-nearest neighbors around each
        georeferenced value that are subsequently kriged.
    wr_indices : |np.ndarray[int]|
        Indices of within-radius (WR) (i.e. < ``k_max``) points.
    oos_indices : |np.ndarray[np.number]|
        Template array based on the size of the data input and ``k_min`` that will contain indices
        where extrapolation is required where there are fewer than ``k_min`` nearest neighbors.
    oos_weights : |np.ndarray[float]|
        Weights applied to extraplolated values.

    Returns
    -------
    Tuple[|np.ndarray[np.number]|, |np.ndarray[np.number]|, |np.ndarray[np.number]|]
        A tuple with updated values for ``wr_indices``, ``oos_indices``, and ``oos_weights`` via a
        search strategy that applies unconstrained and uniform extrapolation to out-of-sample (OOS)
        points.
    """

    # Index for areas with some valid points but fewer than k_min
    partial_indices = sparse_radii[valid_distances[sparse_radii] > 0]

    # Create boolean mask for within-range/sample points
    wr_mask = local_points[partial_indices, :k_max] < search_radius

    # Create temporary matrix for within-range samples
    wr_tmp = wr_indices[partial_indices].copy()
    # ---- Update temporary matrix by applying `wr_mask` for wr points
    wr_tmp[~wr_mask] = np.nan

    # Create temporary matrix for oos samples
    oos_tmp = wr_indices[partial_indices].copy()
    # ---- Update temporary matrix by applying `wr_mask` for oos points
    oos_tmp[wr_mask] = np.nan

    # Assign the OOS values to `oos_indices`
    oos_indices[partial_indices] = np.sort(oos_tmp[:, :k_min])

    # Apply the mask to the remaining `wr_indices` values
    wr_indices[partial_indices] = np.sort(wr_tmp[:, :k_max])

    # Get areas with no valid points within the search radius
    full_extrap_indices = sparse_radii[valid_distances[sparse_radii] == 0]
    if len(full_extrap_indices) > 0:
        # ---- Use all `k_min`-nearest neighbors for extrapolation
        oos_indices[full_extrap_indices] = wr_indices[full_extrap_indices, :k_min]
        wr_indices[full_extrap_indices] = np.nan
        # ---- Compute the OOS kriging weights
        oos_mean = np.apply_along_axis(np.nanmean, 1, local_points[full_extrap_indices, :k_min])
        # ---- Exponentiate the OOS mean
        oos_exp = np.exp(-oos_mean / search_radius)
        # ---- Update the OOS kriging weights
        oos_weights[full_extrap_indices] = oos_exp

    # Return Tuple
    return wr_indices, oos_indices, oos_weights
```
:::
:::{tab-item} Western boundary
:sync: tab2

```python
def western_boundary_search_strategy(
    kriging_mesh: pd.DataFrame,
    western_extent: pd.DataFrame,
    coordinate_names: Tuple[str, str],
    sparse_radii: np.ndarray[int],
    valid_distances: np.ndarray[int],
    local_points: np.ndarray[float],
    distance_matrix_masked: np.ndarray[float],
    nearby_indices: np.ndarray[int],
    k_min: int,
    k_max: int,
    search_radius: float,
    wr_indices: np.ndarray[int],
    oos_indices: np.ndarray[np.number],
    oos_weights: np.ndarray[float],
    **kwargs,
) -> Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]:
    """
    Search strategy that applies western boundary constraints for transect-based surveys

    Parameters
    ----------
    kriging_mesh : pd.DataFrame
        Kriging mesh used for interpolated data values via geostatistics.
    western_extent : pd.DataFrame
        DataFrame with the western-most extent of each transect line used for re-weighting the
        out-of-sample/extrapolated kriged values.
    coordinate_names : Tuple[str, str]
        Names of the coordinate columns when using DataFrames. Expected format: (x_col, y_col).
    sparse_radii : np.ndarray[int]
        Indices where there are fewer than `k_min` nearest neighbors.
    valid_distances : np.ndarray[int]
        The number of masked distance matrix values where extrapolation is required.
    local_points : np.ndarray[float]
        An array with the sorted distances (from nearest to furthest) relative to each point.
    distance_matrix_masked : np.ndarray[float]
        An array with the search-radius-masked nearest neighbor distances.
    nearby_indices : np.ndarray[int]
        Indices of points that require extrapolation.
    k_min : int
        The minimum number of nearest neighbors required for including values for kriging within
        the search radius.
    k_max : int
        The maximum number of nearest neighbors required for including values for kriging detected
        within the search radius.
    search_radius : float
        The adaptive search radius that identifies the *k*-nearest neighbors around each
        georeferenced value that are subsequently kriged.
    wr_indices : np.ndarray[int]
        Indices of within-radius (WR) (i.e. < `k_max`) points.
    oos_indices : np.ndarray[np.number]
        Template array based on the size of the data input and `k_min` that will contain indices
        where extrapolation is required where there are fewer than `k_min` nearest neighbors.
    oos_weights : np.ndarray[float]
        Weights applied to extraplolated values.

    Returns
    -------
    Tuple[np.ndarray[np.number], np.ndarray[np.number], np.ndarray[np.number]]
        A tuple with updated values for `wr_indices`, `oos_indices`, and `oos_weights` via a
        search strategy that uses an extrapolation re-weighting based on transect extents.
    """

    # Parse ordered coordinate names
    x_name, y_name = coordinate_names

    # Index the mesh grid coordinates for bounding the search radius expansion/extrapolation
    # ---- y-coordinates with array transformation to access matrix operations
    mesh_y = kriging_mesh[y_name].to_numpy()[sparse_radii].reshape(-1, 1)
    # ---- x-coordinates
    mesh_x = kriging_mesh[x_name].to_numpy()[sparse_radii]

    # Calculate the mesh distance from the western boundary of the survey transects
    # ---- Find closest point
    mesh_western_distance = np.abs(mesh_y - western_extent["y"].to_numpy()).argmin(axis=1)
    # ---- Calculate the western limits (x-axis)
    western_limit = western_extent.iloc[np.ravel(mesh_western_distance)][x_name]
    # ---- Compute bounding threshold (for tapered extrapolation function)
    western_threshold = western_limit - search_radius
    # ---- Create a thresholded mask for lazy operations
    western_limit_mask = mesh_x < western_threshold

    # Adjust values that don't fall outside the western extent
    if np.any(~western_limit_mask):
        # ---- Grab all values that don't fall outside the western extent
        soft_extrapolation_index = sparse_radii[~western_limit_mask]
        # ---- Find the local points where there are at least some valid points
        partial_indices = soft_extrapolation_index[valid_distances[soft_extrapolation_index] > 0]
        # ---- Update the current values in `wr_indices`
        if len(partial_indices) > 0:
            # -------- Create boolean mask for within-range/sample points
            wr_mask = local_points[partial_indices, :k_max] < search_radius
            # -------- Create temporary matrix for within-range samples
            wr_tmp = wr_indices[partial_indices].copy()
            # -------- Create temporary matrix for oos samples
            oos_tmp = wr_indices[partial_indices].copy()
            # -------- Update temporary matrix by applying `wr_mask` for wr points
            wr_tmp[~wr_mask] = np.nan
            # -------- Update temporary matrix by applying `wr_mask` for oos points
            oos_tmp[wr_mask] = np.nan
            # -------- Assign the OOS values to `oos_indices`
            oos_indices[partial_indices] = np.sort(oos_tmp[:, :k_min])
            # -------- Apply the mask to the remaining `wr_indices` values
            wr_indices[partial_indices] = np.sort(wr_tmp[:, :k_max])

        # ---- Find the local points where there are no valid points within the search radius
        full_extrap_indices = soft_extrapolation_index[
            valid_distances[soft_extrapolation_index] == 0
        ]
        if len(full_extrap_indices) > 0:
            # -------- Update `oos_indices`
            oos_indices[full_extrap_indices] = wr_indices[full_extrap_indices, :k_min]
            # -------- Update `wr_indices`
            wr_indices[full_extrap_indices] = np.nan

    # Taper function for extrapolating values outside the search radius
    if np.any(western_limit_mask):
        # ---- Index these values
        extrapolation_index = sparse_radii[western_limit_mask]
        # ---- Compute the OOS kriging weights
        oos_mean = np.apply_along_axis(np.nanmean, 1, local_points[extrapolation_index, :k_min])
        # ---- Exponentiate the OOS mean
        oos_exp = np.exp(-oos_mean / search_radius)
        # ---- Update the OOS weights
        oos_weights[extrapolation_index] = oos_exp
        # ---- Get the outside indices that correspond to this tapered extrapolation
        sparse_extrapolation_index = nearby_indices[western_limit_mask].astype(float)
        # ---- Apply indices as a mask to the NaN-masked distance matrix
        extrapolated_distance = np.take_along_axis(
            distance_matrix_masked[sparse_radii][western_limit_mask],
            sparse_extrapolation_index.astype(int),
            axis=1,
        )
        # ---- Create NaN mask
        extrapolated_nan_mask = ~np.isnan(extrapolated_distance)
        # -------- Apply mask to indices
        sparse_extrapolation_index_nan = sparse_extrapolation_index.copy()
        sparse_extrapolation_index_nan[extrapolated_nan_mask] = np.nan
        # -------- Update `out_of_sample_indices` matrix
        oos_indices[extrapolation_index] = np.sort(sparse_extrapolation_index_nan)
        # ---- Get inside indices that apply to these points
        # -------- Create NaN mask for within-sample values
        interpolated_nan_mask = np.isnan(extrapolated_distance)
        # -------- Apply mask to indices
        sparse_interpolation_index_nan = sparse_extrapolation_index.copy()
        sparse_interpolation_index_nan[interpolated_nan_mask] = np.nan
        # -------- Pad NaN to match `within_sample_indices` matrix
        sparse_interpolation_pad = np.pad(
            sparse_interpolation_index_nan,
            [(0, 0), (0, k_max - k_min)],
            mode="constant",
            constant_values=np.nan,
        )
        # -------- Updated `within_sample_indices` matrix
        wr_indices[extrapolation_index] = np.sort(sparse_interpolation_pad)

    # Return Tuple
    return wr_indices, oos_indices, oos_weights
```
:::
::::