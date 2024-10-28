from pydantic import BaseModel, ConfigDict, Field, ValidationError, field_validator
from typing import Dict, Literal, Optional, Tuple
import cartopy.feature as cfeature
import numpy as np
import re

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
    
    figure_width: float = Field(default=5.5, gt=0.0, allow_inf_nan=False)
    axis_limits: Optional[Dict[str, Tuple[float, float]]] = Field(default=None)
    colormap: Optional[str] = Field(default=None)
    data_range: Optional[Tuple[float, float]] = Field(default=None)
    log_base: Optional[float] = Field(default=None, gt=0.0, allow_inf_nan=False)
    model_config = ConfigDict(extra="allow")
       
    @field_validator("data_range", mode="before")
    def validate_data_range(cls, v):
        
        # Check whether `v` is None or is a valid tuple
        if v is None:
            return v
        
        # Validate the expected 2-tuple format
        if not isinstance(v, tuple) or len(v) != 2:
            raise ValueError(
                "Input for `data_range` must comprise a tuple with exactly exactly 2 float values."
                )
            
        # Try coercing to a float
        try:
            v_tpl = tuple(float(i) for i in v)
        except Exception as e:
            raise e
            
        # Validate that values are floats and are finite
        if not all(isinstance(i, float) and np.isfinite(i) for i in v_tpl):
            raise ValueError(
                "Input values for `data_range` must all be real numbers."
                )
            
        # Return coerced value
        return v_tpl

    @field_validator("axis_limits", mode="before")
    def validate_axis_limits(cls, v):
        
        # Check whether `v` is None or is a valid tuple
        if v is None:
            return v
        
        # Validate the expected dictionary format
        if not isinstance(v, dict):
            raise ValueError(
                "Input for `axis_limits` must be a dictionary."
            )
        
        # Validate all items within the dictionary
        for key, value in v.items():
            # ---- Validate the expected 2-tuple format
            if not isinstance(value, tuple) or len(value) != 2:
                raise ValueError(
                    f"Input for `axis_limits[{key}]` must comprise a tuple with exactly 2 float "
                    f"values."
                    )
            # ---- Try coercing to a float
            try:
                v_tpl = tuple(float(i) for i in value)
            except Exception as e:
                raise e
            # ---- Validate that values are both floats and finite
            if not all(isinstance(i, float) and np.isfinite(i) for i in v_tpl):
                raise ValueError(
                    f"Input for `axis_limits[{key}]` must all be real numbers."
                )             

    # Factory method
    @classmethod
    def create(cls, **kwargs):
        """
        Factory creation method
        """

        # Scrutinize `kind`
        cls.judge(**kwargs)
        
        # Find the correct sub-class plotting model
        return (
            cls._SUBCLASS_FACTORY()[kwargs["kind"]]
            .judge(**kwargs).model_dump(exclude_none=False)
        )
    
    @classmethod
    def _SUBCLASS_FACTORY(cls):
        
        return {
            "age_length_distribution": BiologicalHeatmapPlot,
            "mesh": MeshPlot,
            "transect": TransectPlot,
        }
          
class GeoConfigPlot(BaseModel, arbitrary_types_allowed=True):
    """
    Geospatial parameters

    Parameters
    ----------
    init: str
        EPSG projection code.
    """

    init: str = Field(default="epsg:4326", description="EPSG projection code.")
    coastline: cfeature.NaturalEarthFeature = Field(
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
    
    geo_config: GeoConfigPlot = Field(default_factory=lambda: GeoConfigPlot.create())    

class TransectPlot(SpatialPlot):
    
    variable: Literal["abundance", "abundance_female", "abundance_male", "biomass", 
                      "biomass_density", "biomass_density_female", "biomass_density_male", 
                      "biomass_female", "biomass_male", "nasc", "number_density", 
                      "number_density_female", "number_density_male"]
    
class MeshPlot(SpatialPlot):
    
    variable: Literal["biomass", "kriged_mean", "kriged_variance", "sample_cv", "sample_variance"]

class BiologicalHeatmapPlot(PlotModel):
    """
    Base Pydantic model for biological heatmap plots
    """
    
    variable: Literal["abundance", "biomass"]
    sex: Literal["all", "female", "male"] = Field(default="all")
    grid: bool = Field(default=True)