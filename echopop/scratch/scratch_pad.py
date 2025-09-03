data = df_nasc_no_age1
import param

class CompleteVariogramGUI(param.Parameterized):
    """Complete interactive GUI for variogram analysis with base parameter controls"""
    
    def __init__(self, data=None, **params):
        super().__init__(**params)
        self.data = data
        self.vgm = None  # Will be created based on GUI parameters
        self.empirical_results = {}
        self.model_results = {}
        self.optimization_results = {}
        
self = CompleteVariogramGUI(data)
