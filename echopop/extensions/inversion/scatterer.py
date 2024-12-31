"""
Scatterer class used for parameterizing and computing acoustic scattering models
"""

from typing import Any, Dict, Literal

class Scatterer:
    
    def __init__(
        self,
        parameters: Dict[str, Any],
        scattering_type: Literal[
            "elastic_shelled",
            "fluid_like",
            "resonant",
        ] = "fluid_like",
        shape: Literal[
            "arbitrary",
            "cylinder",
            "oblate_spheroid",
            "prolate_spheroid",
            "sphere",
            "umbrella"            
        ] = "arbitrary"
    ):
        
        # Initialize attributes
        # ---- Model parameters
        self.parameters = parameters
        # ---- Scattering-type
        self.scattering_type = scattering_type 
        # ---- Shape-type
        self.shape = shape
