"""
This script establishes all global variables used throughout EchoPro.
"""

from typing import TypedDict, Callable


# define the semi-variogram input types
vario_type_dict = {'nlag': int, 'lag_res': float}
vario_param_type = TypedDict('vario_param_type', vario_type_dict)

# define the Kriging parameter input types
krig_type_dict = {'k_max': int, 'k_min': int, 'R': float, 'ratio': float,
                  's_v_params': dict, 's_v_model': Callable}
krig_param_type = TypedDict('krig_param_type', krig_type_dict)
