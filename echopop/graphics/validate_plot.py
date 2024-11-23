import re
from typing import Any, Callable, Dict, Literal, Optional, Union

import cartopy.feature as cfeature
import numpy as np
from cartopy.crs import PlateCarree, Projection
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    SerializeAsAny,
    ValidationError,
    field_validator,
    model_validator,
)


class BasePlotModel(BaseModel):
    """
    Base Pydantic model for scrutinizing plotting parameters
    """

    kind: Literal["age_length_distribution", "mesh", "transect"]

    # Validator method
    @classmethod
    def judge(cls, **kwargs):
        """
        Validator method
        """
        try:
            return cls(**kwargs)
        except ValidationError as e:
            e.__traceback__ = None
            raise e


class PlotModel(BasePlotModel):
    """
    Base Pydantic model for scrutinizing plotting parameters
    """

    axis_limits: Optional[Dict[str, Any]] = Field(default=None)
    log_base: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    model_config = ConfigDict(extra="allow")

    @field_validator("axis_limits", mode="before")
    def validate_axis_limits(cls, v):

        # Check whether `v` is None or is a valid tuple
        if v is None:
            return v

        # Validate the expected dictionary format
        if not isinstance(v, dict):
            raise ValueError("Input for `axis_limits` must be a dictionary.")

        # Validate all items within the dictionary
        for key, value in v.items():
            # ---- Validate the expected 2-value format
            if set(["xmin", "xmax"]).issubset(value):
                test_values = value
            elif set(["ymin", "ymax"]).issubset(value):
                test_values = value
            elif set(["left", "right"]).issubset(value):
                test_values = value
            else:
                raise ValueError(
                    f"Input for `axis_limits['{key}']` must comprise either two values represented "
                    f"represented by '{key}min'/'{key}max' or 'left'/'right'."
                )
            # ---- Check that both values are floats
            if not all([isinstance(val, (float, int)) for _, val in test_values.items()]):
                raise ValueError(
                    f"Input for `axis_limits['{key}']` must comprise float values for limits."
                )

        # Return
        return v

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """

        # Scrutinize `kind`
        cls.judge(**kwargs)

        # Find the correct sub-class plotting model (plys validation)
        # ---- Get the specific plot parameters based on `kind`
        kind_params = (
            cls._SUBCLASS_FACTORY()[kwargs["kind"]].judge(**kwargs).model_dump(exclude_none=False)
        )
        # ---- Get the specific plot parameters based on `plot_type` and `variable`
        plot_params = ReferenceParams.create(**{**kwargs, **kind_params})

        # ---- Combine and return
        return {**plot_params, **kind_params}

    @classmethod
    def _SUBCLASS_FACTORY(cls):

        return {
            "age_length_distribution": BiologicalHeatmapPlot,
            "mesh": MeshPlot,
            "transect": TransectPlot,
        }


# --------------------------------------------------------------------------------------------------
# Plot-kind validators
class GeoConfigPlot(BaseModel, arbitrary_types_allowed=True):
    """
    Geospatial parameters

    Parameters
    ----------
    init: str
        EPSG projection code.
    coastline: cartopy.feature.NaturalEarthFeature
        A land or coastline feature used for geospatial plotting.
    """

    coastline: cfeature.Feature = Field(
        default_factory=lambda: cfeature.NaturalEarthFeature(
            category="physical",
            name="land",
            scale="10m",
            facecolor="#C5C9C7",
            edgecolor="#929591",
            linewidth=0.5,
            zorder=1,
        )
    )
    init: str = Field(default="epsg:4326", description="EPSG projection code.")
    plot_projection: Projection = Field(default=PlateCarree())
    model_config = ConfigDict(extra="allow")

    @field_validator("init", mode="before")
    def validate_init(cls, v):
        # ---- Convert to a string if read in as an integer
        if isinstance(v, (int, float)):
            v = str(v)
        # ---- Convert to lowercase
        v = v.lower()
        # ---- Mold the entry into the expected format that includes a preceding 'epsg:'
        if not v.startswith("epsg"):
            v = "epsg:" + v
        # ---- Ensure that the colon is present
        if ":" not in v:
            v = "epsg:" + v.split("epsg")[1]
        # ---- Evaluate whether the pre-validator succeeded in finding an acceptable format
        if not re.match(r"^epsg:\d+$", v):
            raise ValueError(
                f"Echopop cannot parse the defined EPSG code ('{v}'). EPSG codes most be formatted "
                f"with strings beginning with 'epsg:' followed by the integer number code (e.g. "
                f"'epsg:4326')."
            )
        # ---- Return the pre-validated entry
        return v

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """

        return cls(**kwargs).model_dump(exclude_none=False)


class SpatialPlot(PlotModel):
    """
    Base Pydantic model for spatial plots
    """

    geo_config: SerializeAsAny[GeoConfigPlot] = Field(
        default_factory=lambda: GeoConfigPlot.create()
    )


class TransectPlot(SpatialPlot):

    plot_type: Literal["scatter"] = Field(default="scatter")
    variable: Literal[
        "abundance",
        "abundance_female",
        "abundance_male",
        "biomass",
        "biomass_density",
        "biomass_density_female",
        "biomass_density_male",
        "biomass_female",
        "biomass_male",
        "nasc",
        "number_density",
        "number_density_female",
        "number_density_male",
    ] = Field(default="biomass")

    @field_validator("plot_type", mode="before")
    def validate_plot_type(cls, v):
        if v is None:
            return "scatter"
        else:
            return v


class MeshPlot(SpatialPlot):

    plot_type: Optional[Literal["hexbin", "pcolormesh", "scatter"]] = Field(default="hexbin")
    variable: Literal[
        "biomass", "biomass_density", "kriged_variance", "kriged_cv", "local_variance"
    ] = Field(default="biomass")

    @field_validator("plot_type", mode="before")
    def validate_plot_type(cls, v):
        if v is None:
            return "hexbin"
        else:
            return v


class BiologicalHeatmapPlot(PlotModel):
    """
    Base Pydantic model for biological heatmap plots
    """

    grid_heatmap: bool = Field(default=True)
    plot_type: Literal["heatmap"] = Field(default="heatmap")
    sex: Literal["all", "female", "male"] = Field(default="all")
    variable: Literal["abundance", "biomass"] = Field(default="biomass")

    @field_validator("plot_type", mode="before")
    def validate_plot_type(cls, v):
        if v is None:
            return "heatmap"
        else:
            return v


# --------------------------------------------------------------------------------------------------
# Plot-type validators
class BaseTypeVar(BaseModel):
    model_config = ConfigDict(extra="allow")

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        try:
            return cls(**kwargs).model_dump(exclude_none=True)
        except ValidationError as e:
            e.__traceback__ = None
            raise e


class HeatmapVar(BaseTypeVar):
    aspect: Union[Literal["auto", "equal"], float] = Field(default="auto")


class HexbinVar(BaseTypeVar):
    linewidth: float = Field(
        default=0.05, description="`matplotlib.pyplot.hexbin` argument for line widths."
    )
    reduce_C_function: Optional[Callable] = Field(
        default=None,
        description="Function to use for reducing the \
                                                      C-dimension.",
    )
    variable: str

    @model_validator(mode="before")
    def validate_reduce_C_function(cls, v):

        # Get `variable`
        variable = v.get("variable")

        # Get current `reduce_C_function`
        if not v.get("reduce_C_function", None):
            if variable == "biomass":
                v["reduce_C_function"] = np.sum
            else:
                v["reduce_C_function"] = np.mean

        # Return
        return v


class ScatterVar(BaseTypeVar):
    edgecolors: Optional[str] = Field(
        default=None,
        description="`matplotlib.pyplot.scatter` argument for \
                                         edgecolor",
    )
    linewidths: Optional[float] = Field(
        default=None,
        description="`matplotlib.pyplot.scatter` argument for \
                                           linewidth",
    )
    marker: Optional[str] = Field(
        default=None,
        description="`matplotlib.pyplot.scatter` argument for marker \
                                      shape.",
    )
    s: Optional[Union[Callable, float]] = Field(
        default=None, description="`matplotlib.pyplot.scatter` argument for marker size."
    )
    vmax: Optional[float] = Field(default=None)
    vmin: Optional[float] = Field(default=None)

    @model_validator(mode="before")
    def set_s_default(cls, v):

        # Only apply if no value is supplied for `s` and `kind` is "transect"
        if v.get("s", None) is None and v.get("kind") == "transect":
            # ---- Create wrapper
            def scale_sizes_wrapper(x, vmin: float, vmax: float):
                return scale_sizes(x, vmin, vmax, 1, 75)

            # ---- Create scaling Callable
            v["s"] = scale_sizes_wrapper
        elif v.get("s", None) is None:
            # ---- Otherwise, set static
            v["s"] = 2.0

        # Only apply if no value is supplied for `marker` and `kind` is "mesh"
        if v.get("kind") == "mesh":
            # ---- Create scaling Callable
            v["marker"] = v.get("marker", "s")

        # Only apply if no values are supplied relevant for `kind="transect"`
        if v.get("kind") == "transect":
            # ---- `edgecolor`
            v["edgecolors"] = v.get("edgecolors", "black")
            # ---- `linewidth`
            v["linewidths"] = v.get("linewidths", 0.2)

        # Return value
        return v


class PcolormeshVar(BaseTypeVar):
    add_colorbar: bool = Field(default=False)


# --------------------------------------------------------------------------------------------------
# Plotting parameter validators
class BaseVarRef(BaseModel):
    vmin: float = Field(default=0.0, description="Minimum value for the colormap.")
    model_config = ConfigDict(extra="allow")

    _is_colorbar_label_present: bool = PrivateAttr(default=False, init=False)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        if "colorbar_label" in kwargs:
            self._is_colorbar_label_present = True

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """
        try:
            return cls(**kwargs).model_dump(exclude_none=True)
        except ValidationError as e:
            e.__traceback__ = None
            raise e

    @model_validator(mode="after")
    def adjust_colorbar_label(self):

        if not self._is_colorbar_label_present:

            # Get the variable
            variable = self.variable

            # Get kind of plot
            kind = self.kind

            # Add prepend label string, if required
            if kind in ["age_length_distribution", "transect"]:
                if "_male" in variable:
                    label_prepend = "Male "
                elif "_female" in variable:
                    label_prepend = "Female "
                elif variable in ["abundance", "biomass", "biomass_density", "number_density"]:
                    label_prepend = ""
                else:
                    label_prepend = None
            elif kind == "mesh":
                if variable not in ["kriged_cv", "kriged_variance", "local_variance"]:
                    label_prepend = "Kriged "
                else:
                    label_prepend = None

            # Prepend the determined label to the default value
            if label_prepend:
                self.colorbar_label = (label_prepend + self.colorbar_label).capitalize()

        return self


class AbundanceVar(BaseVarRef, arbitrary_types_allowed=True):
    cmap: str = Field(default="viridis")
    colorbar_label: str = Field(default="Abundance\n# animals")
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class BiomassVar(BaseVarRef):
    cmap: str = Field(default="inferno")
    colorbar_label: str = Field(default="Biomass\n$\\mathregular{kg}$")
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class BiomassDensityVar(BaseVarRef):
    cmap: str = Field(default="plasma")
    colorbar_label: str = Field(default="Biomass density\n$\\mathregular{kg~nmi^{-2}}$")
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class KrigedCVVar(BaseVarRef):
    cmap: str = Field(default="magma")
    colorbar_label: str = Field(default="Kriged $CV$")
    vmax: Union[float, Callable] = Field(default=lambda v: np.ceil(v.max() / 0.1) * 0.1)


class KrigedVarianceVar(BaseVarRef):
    cmap: str = Field(default="hot")
    colorbar_label: str = Field(
        default="Kriged biomass density variance\n$\\mathregular{(kg~nmi^{-2}})^{2}$"
    )
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class LocalVarianceVar(BaseVarRef):
    cmap: str = Field(default="cividis")
    colorbar_label: str = Field(
        default="Local biomass density variance\n$\\mathregular{(kg~nmi^{-2}})^{2}$"
    )
    reduce_C_function: Callable = Field(default=np.mean)
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class NASCVar(BaseVarRef):
    cmap: str = Field(default="cividis")
    colorbar_label: str = Field(default="NASC\n$\\mathregular{m^{2}~nmi^{-2}}$")
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class NumberDensityVar(BaseVarRef):
    cmap: str = Field(default="magma")
    colorbar_label: str = Field(default="Number density\n$\\mathregular{animals~nmi^{-2}}$")
    vmax: Union[float, Callable] = Field(default=lambda v: 10 ** np.round(np.log10(v.max())))


class ReferenceParams(BaseVarRef):
    cmap: Optional[str] = Field(default=None, description="Colormap name.")
    colorbar_label: Optional[str] = Field(default=None, description="Colorbar label.")
    plot_type: str
    variable: str
    vmax: Optional[Union[float, Callable]] = Field(
        default=None, description="Maximum value for the colormap."
    )
    vmin: Optional[float] = Field(default=None, description="Minimum value for the colormap.")

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """

        # Get plot-based parameters (plus validation)
        type_params = cls._REFERENCE_TYPE_FACTORY(kwargs["plot_type"]).create(**kwargs)

        # Get variable-based parameters (plus validation)
        var_params = cls._REFERENCE_VAR_FACTORY(kwargs["variable"]).create(**kwargs)

        # Combine and return
        return {**type_params, **var_params}

    @classmethod
    def _REFERENCE_VAR_FACTORY(cls, variable: str):

        if "abundance" in variable:
            return AbundanceVar
        elif "biomass_density" in variable:
            return BiomassDensityVar
        elif "biomass" in variable:
            return BiomassVar
        elif "nasc" in variable:
            return NASCVar
        elif "number_density" in variable:
            return NumberDensityVar
        elif "local_variance" in variable:
            return LocalVarianceVar
        elif "kriged_variance" in variable:
            return KrigedVarianceVar
        elif "kriged_cv" in variable:
            return KrigedCVVar

    @classmethod
    def _REFERENCE_TYPE_FACTORY(cls, plot_type: str):

        if plot_type == "hexbin":
            return HexbinVar
        elif plot_type == "scatter":
            return ScatterVar
        elif plot_type == "pcolormesh":
            return PcolormeshVar
        elif plot_type == "heatmap":
            return HeatmapVar


# --------------------------------------------------------------------------------------------------
# Utility function


def scale_sizes(values, min_value, max_value, min_size=25, max_size=250):
    """
    Scale point size
    """
    # Censor values if needed
    sizes = values.copy()
    sizes.loc[sizes < min_value] = min_value
    sizes.loc[sizes > max_value] = max_value

    return ((sizes - min_value) / (max_value - min_value)) * (max_size - min_size) + min_size
