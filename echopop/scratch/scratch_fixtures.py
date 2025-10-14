import lmfit
import numpy as np
from lmfit import minimize, Minimizer, fit_report

class MyParameter(lmfit.Parameter):
    """A custom Parameter class that stores original min and max for unscaling."""
    def __init__(self, *args, original_min=None, original_max=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.original_min = original_min
        self.original_max = original_max

class MyParameters(lmfit.Parameters):
    """A custom Parameters class that normalizes inputs and stores MyParameter instances."""
    def add(self, name, value=None, vary=True, min=None, max=None, expr=None, brute_step=None):
        if isinstance(name, MyParameter):
            self[name.name] = name
        else:
            # Check for min and max to perform normalization
            if min is not None and max is not None:
                original_min = min
                original_max = max
                normalized_value = (value - original_min) / (original_max - original_min)
                self[name] = MyParameter(name=name, value=normalized_value, vary=vary,
                                         min=0.0, max=1.0, expr=expr,
                                         brute_step=brute_step,
                                         original_min=original_min,
                                         original_max=original_max)
            else:
                # If no min/max provided, add a MyParameter with None for original bounds
                self[name] = MyParameter(name=name, value=value, vary=vary,
                                         min=min, max=max, expr=expr,
                                         brute_step=brute_step,
                                         original_min=None,
                                         original_max=None)
                
def objective(my_parameters, data):
    """
    Objective function that receives normalized parameters and un-normalizes them
    for the model calculation.
    """
    a_param = my_parameters['a']
    b_param = my_parameters['b']

    # Perform the correct unscaling operation if bounds exist
    if a_param.original_min is not None and a_param.original_max is not None:
        a_unscaled = a_param.value * (a_param.original_max - a_param.original_min) + a_param.original_min
    else:
        a_unscaled = a_param.value # Use the value directly if not normalized

    if b_param.original_min is not None and b_param.original_max is not None:
        b_unscaled = b_param.value * (b_param.original_max - b_param.original_min) + b_param.original_min
    else:
        b_unscaled = b_param.value # Use the value directly if not normalized

    # Use the unscaled values for the model calculation
    return a_unscaled * data['x'] + b_unscaled - data['y']

# Generate some data
x_data = np.linspace(0, 10, 100)
y_data = 2.5 * x_data + 1.2 + np.random.normal(0, 0.5, size=x_data.shape)
fit_data = {'x': x_data, 'y': y_data}

# Now, we provide the initial values and bounds in their original units.
# The `add` method will handle the normalization internally.
params = MyParameters()
params.add('a', value=2.0, min=0.0, max=5.0)
params.add('b', value=1.0, min=0.0, max=2.0)

# Create a Minimizer and run the fit.
mini = Minimizer(objective, params, fcn_kws=dict(data=fit_data))
result = mini.minimize(method='least_squares')

# Print the report
print(fit_report(result))