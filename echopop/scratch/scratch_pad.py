VARIOGRAM_MODEL_PARAMETER_MAP = {
    "Bessel-Exponential": ["bessel", "exponential"],
    "Bessel-Gaussian": ["bessel", "gaussian"],
    "Cosine-Exponential": ["cosine", "exponential"],
    "Cosine-Gaussian": ["cosine", "gaussian"],
    "Cubic": "cubic",
    "Exponential": "exponential",
    "Exponential-Linear": ["exponential", "linear"],
    "Gaussian": "gaussian",
    "Gaussian-Linear": ["gaussian", "linear"],
    "J-Bessel": "jbessel",
    "K-Bessel": "kbessel",
    "Linear": "linear",
    "Matérn": "matern",
    "Nugget": "nugget",
    "Pentaspherical": "pentaspherical",
    "Power law": "power",
    "Rational quadratic": "quadratic",
    "Cardinal sine": "sinc",
    "Spherical": "spherical",
}

DEFAULT_VARIOGRAM_PARAMETERS = {
    "correlation_range": {
        "name": r"Correlation range ($a$)" ,
        "min": 1e-10, "value": 1.0, "max": 99999, "vary": False, "step": 1.0
    },
    "decay_power": {
        "name": r"Decay power ($\omega$)",
        "min": 1e-10, "value": 1.0, "max": 2.0, "vary": False, "step": 0.1
    },
    "enhance_semivariance": {
        "name": "Enhance semivariance", 
        "value": True
    },
    "hole_effect_range": {
        "name": r"Hole effect range ($a_h$)",
        "min": 0.0, "value": 0.0, "max": 99999, "vary": False, "step": 1.0
    },
    "sill": {
        "name": r"Sill ($C$)",
        "min": 1e-10, "value": 1.0, "max": 99999, "vary": False, "step": 1.0
    },
    "nugget": {
        "name": r"Nugget ($C_0$)",
        "min": 0.0, "value": 0.0, "max": 99999, "vary": False, "step": 1.0
    },
    "smoothness_parameter": {
        "name": r"Matérn shape parameter ($\nu$)",
        "min": 0.0, "value": 0.5, "max": 10.0, "vary": False, "step": 0.1
    },
    "shape_parameter": {
        "name": r"Scale ($\beta$)",
        "min": 1e-10, "value": 1.0, "max": 100.0, "vary": False, "step": 1.0
    },
    "power_exponent": {
        "name": r"Power ($\omega$)",
        "min": 1e-10, "value": 1.0, "max": 2.0 - 1e-10, "vary": False, "step": 0.1
    }
}

VARIABLE_MAP = {
    "Number density (animals nmi^-2)": "number_density", 
    "Biomass density (kg nmi^-2)": "biomass_density", 
    "NASC (m^2 nmi^-2)": "nasc",
}


class CompleteVariogramGUI(param.Parameterized):
    """Complete interactive GUI for variogram analysis with base parameter controls"""
    
    def __init__(self, data=None, **params):
        super().__init__(**params)
        self.data = data
        self.vgm = None  # Will be created based on GUI parameters
        self.empirical_results = {}
        self.model_results = {}
        self.optimization_results = {}

import holoviews as hv
hv.extension('bokeh')
from bokeh.models import HoverTool
from IPython.display import display_html
self = CompleteVariogramGUI(data)

self.vgm = vgm


####
lags, gamma, lag_counts = (self.vgm.lags, self.vgm.gamma, self.vgm.lag_counts)

data = pd.DataFrame({
    "distance": lags, 
    "semivariance": gamma, 
    "lag_count": lag_counts, 
    "lag_number": range(1, len(lags) + 1)
}) 

# Scale point sizes based on lag counts
if len(lag_counts) > 1 and lag_counts.max() > lag_counts.min():
    size_scale = lambda x: ((x - x.min()) / (x.max() - x.min()) + 1) * 15
    data['point_size'] = size_scale(lag_counts)
else:
    data['point_size'] = 15

scatter = hv.Scatter(
    data, 
    kdims=['distance'], 
    vdims=['semivariance', 'lag_count', 'lag_number', 'point_size']
).opts(
    size='point_size',
    color='blue',
    line_color='black',
    line_width=1,
    width=700,
    height=400,
    xlabel='Lag distance [h]',
    ylabel='Semivariance [γ]',
    title='Empirical Variogram',
    show_grid=True,
    # Explicitly specify tools to prevent duplication and include reset
    tools=['hover', 'pan', 'wheel_zoom', 'box_zoom', 'reset'],
    active_tools=['pan']
)

# Custom hover tool configuration
hover = HoverTool(
    tooltips=[
        ('Lag #', '@lag_number'),
        ('Distance', '@distance{0.000}'),
        ('Semivariance', '@semivariance{0.000}'),
        ('Count', '@lag_count')
    ]
)

hv.output(scatter)
