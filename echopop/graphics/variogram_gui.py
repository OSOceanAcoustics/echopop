import ast
import traceback
from typing import Literal, Optional

import holoviews as hv
import ipywidgets as ipw
import markdown
import pandas as pd
from bokeh.models import HoverTool
from IPython.display import clear_output, display
from lmfit import Parameters

from ..geostatistics import Variogram, compute_variogram, get_variogram_arguments

# ==================================================================================================
# Constants, mappings, and helper functions API
# ---------------------------------------------
# Vertical and horizontal spacers
VSPACER = ipw.Box(layout=ipw.Layout(height="5px"))
HSPACER = ipw.Box(layout=ipw.Layout(width="10px"))

# Variogram variable mapping
VARIABLE_MAP = {
    "Number density (animals nmi⁻²)": "number_density",
    "Biomass density (kg nmi⁻²)": "biomass_density",
    "NASC (m² nmi⁻²)": "nasc",
}

# Variogram model mapping
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

# Default parameter values mapping
DEFAULT_VARIOGRAM_PARAMETERS = {
    "correlation_range": {
        "name": "Correlation range (a)",
        "short_name": "a",
        "min": 1e-10,
        "value": 1.0,
        "max": 99999,
        "vary": True,
        "step": 1.0,
    },
    "decay_power": {
        "name": "Decay power (\u03b1)",
        "short_name": "\u03b1",
        "min": 1e-10,
        "value": 1.0,
        "max": 2.0,
        "vary": True,
        "step": 0.1,
    },
    "enhance_semivariance": {"name": "Enhance semivariance", "value": True},
    "hole_effect_range": {
        "name": "Hole effect range (a\u2095)",
        "short_name": "a\u2095",
        "min": 0.0,
        "value": 0.0,
        "max": 99999,
        "vary": True,
        "step": 1.0,
    },
    "sill": {
        "name": "Sill (C)",
        "short_name": "C",
        "min": 1e-10,
        "value": 1.0,
        "max": 99999,
        "vary": True,
        "step": 1.0,
    },
    "nugget": {
        "name": "Nugget (C\u2080)",
        "short_name": "C\u2080",
        "min": 0.0,
        "value": 0.0,
        "max": 99999,
        "vary": True,
        "step": 1.0,
    },
    "smoothness_parameter": {
        "name": "Matérn shape parameter (\u03bd)",
        "short_name": "\u03bd",
        "min": 0.0,
        "value": 0.5,
        "max": 10.0,
        "vary": True,
        "step": 0.1,
    },
    "shape_parameter": {
        "name": "Scale (\u03b2)",
        "short_name": "\u03b2",
        "min": 1e-10,
        "value": 1.0,
        "max": 100.0,
        "vary": True,
        "step": 1.0,
    },
    "power_exponent": {
        "name": "Power (\u03c9)",
        "short_name": "\u03c9",
        "min": 1e-10,
        "value": 1.0,
        "max": 2.0 - 1e-10,
        "vary": True,
        "step": 0.1,
    },
}


def md_to_HTML(
    text: str,
    type: Optional[Literal["fail", "status", "success"]] = None,
    ipywidget_obj: Optional[ipw.HTML] = None,
):
    """
    Helper function that converts text and regular expressions from a markdown format to a
    `ipywidgets.HTML` text object. When an `ipywidgets.HTML` object is supplied, then the text is
    stored directly into the 'value' attribute of that object.
    """

    # Get background color
    if type:
        if type == "fail":
            background_color = "lightcoral"
        elif type == "status":
            background_color = "lightblue"
        elif type == "success":
            background_color = "lightgreen"
    else:
        background_color = "white"

    # Stylize text
    styled_text = (
        f"<div style='background-color: {background_color}; padding: 10px; "
        f"border-radius: 5px;'>{markdown.markdown(text)}</div>"
    )

    # Update
    if ipywidget_obj:
        ipywidget_obj.value = markdown.markdown(styled_text)
    else:
        return ipw.HTML(value=markdown.markdown(styled_text))


# ==================================================================================================
# Instructions tab
# ----------------
def instructions_tab(obj):
    instructions_text = f"""
# <u>Variogram Analysis Interactive GUI</u>

<span style='font-size:24px; line-height:1; margin:0; padding:0; display:block;'>
<b>Overview</b></span>

<span style='font-size:18px; line-height:1; margin:0; padding:0; display:block;'>
<b><i>Variogram initialization</i></b></span>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Set the values for the <b>lag resolution</b> and <b>number of lags</b> used for computing the
empirical and theoretical variograms. The initial values correspond to those used for creating
the <span style='color: #4EC9B0; background-color: white; padding: 2px 4px; border-radius: 3px;
font-family: monospace; font-weight: bold;'>VariogramGUI</span> object:
</span>

<ul style='margin: 5px 0; padding-left: 20px; line-height: 1.2;'>
<li style='margin: 2px 0; padding: 0;'><b>Lag resolution:</b>
{obj.initialization_args['lag_resolution']}</li>
<li style='margin: 2px 0; padding: 0;'><b>Number of lags:</b>
{obj.initialization_args['n_lags']}</li>
<li style='margin: 2px 0; padding: 0;'><b>Coordinates:</b>
{obj.initialization_args['coordinates']}</li>
</ul>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
The coordinates correspond to the <i>x</i>- and <i>y</i>-coordinates of the geospatial dataset.
This variable is static and remains constant. Once the variables have been changed to the
desired values, click <span style='background-color: #337ab7; padding: 2px 4px;
border-radius: 3px; font-family: monospace; font-weight: bold; color: white;'>Initialize
variogram</span> and proceed to the next tab: <span style='background-color: lightgrey;
color: black; padding: 2px 4px; border-radius: 3px; font-family: monospace; font-weight: bold;'>
Empirical variogram</span>.
</span>

<span style='font-size:18px; line-height:1; margin:0; padding:0; display:block;'>
<b><i>Empirical variogram</i></b></span>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Select the variable that will be used for the variogram analysis. There are currently three
options available:</span>

<ul style='margin: 5px 0; padding-left: 20px; line-height: 1.2;'>
    {''.join(
        [f"<li style='margin: 2px 0; padding: 0;'><b>{key}</b></li>" for key in VARIABLE_MAP.keys()]
    )}
</ul>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
The azimuth angle filter is optional and will only be applied if the box is checked. This
accounts for directionality in spatial autocorrelation. Valid values are bounded between 0° and
180°. Once the parameters are set to their desired values, click <span
style='background-color: #337ab7; padding: 2px 4px; border-radius: 3px; font-family: monospace;
font-weight: bold; color: white;'>Compute empirical variogram</span>. This will calculate the
semivariance (𝜸) across the defined lag distances (𝒉), which will render an interactive figure
to the right of the control panel. Hovering the mouse cursor over individual points, whose
sizes are scaled to their respective lag counts, will provide lag-specific details. Once this
step is completed, proceed to the next tab: <span style='background-color: lightgrey;
color: black; padding: 2px 4px; border-radius: 3px; font-family: monospace; font-weight: bold;'>
Variogram model</span>.
</span>
</span>

<span style='font-size:18px; line-height:1; margin:0; padding:0; display:block;'>
<b><i>Variogram model</i></b></span>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Select one of the specific variogram models available in the <b>`Variogram model`</b> dropdown
menu. The current models are available:
</span>

<div style='display: flex; margin: 5px 0;'>
    <ul style='margin: 0; padding-left: 20px; padding-right: 5px; line-height: 1.2; width: 50%;'>
    {''.join([
        f"<li style='margin: 2px 0; padding: 0;'><b>{key}</b></li>"
        for key in list(VARIOGRAM_MODEL_PARAMETER_MAP.keys())[
            :len(VARIOGRAM_MODEL_PARAMETER_MAP) // 2
        ]
    ])}
    </ul>
    <ul style='margin: 0; padding-left: 5px; line-height: 1.2; width: 50%;'>
    {''.join([
        f"<li style='margin: 2px 0; padding: 0;'><b>{key}</b></li>"
        for key in list(VARIOGRAM_MODEL_PARAMETER_MAP.keys())[
            len(VARIOGRAM_MODEL_PARAMETER_MAP) // 2:
        ]
    ])}
    </ul>
</div>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Each model selection from the dropdown menu will dynamically populate the window with the model-
specific parameters. For example, selecting <b>`Exponential`</b> will bring up the following
parameters:
</span>

<ul style='margin: 5px 0; padding-left: 20px; line-height: 1.2;'>
    {''.join([
        f"<li style='margin: 2px 0; padding: 0;'><b>{v['name']}</b></li>"
        for k, v in DEFAULT_VARIOGRAM_PARAMETERS.items()
        if k in ['correlation_range', 'nugget', 'sill']
    ])}
</ul>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Alternatively, selecting <b>`Matérn`</b> produces:
</span>

<ul style='margin: 5px 0; padding-left: 20px; line-height: 1.2;'>
    {''.join([
        f"<li style='margin: 2px 0; padding: 0;'><b>{v['name']}</b></li>"
        for k, v in DEFAULT_VARIOGRAM_PARAMETERS.items()
        if k in ['correlation_range', 'nugget', 'sill', 'smoothness_parameter']
    ])}
</ul>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Once the parameters are set to their desired values, click <span
style='background-color: #337ab7; padding: 2px 4px; border-radius: 3px; font-family: monospace;
font-weight: bold; color: white;'>Compute theoretical variogram</span>. This will model the
semivariance and overlay a <span style='background-color:black; color:white;'><b>black line</b>
</span> on top of the empirical variogram results. When a dictionary is provided for <span
style='background-color:white; font-family: monospace'>variogram_parameters</span> when
creating the <span style='color: #4EC9B0; background-color: white; padding: 2px 4px;
border-radius: 3px; font-family: monospace; font-weight: bold;'> VariogramGUI</span> object,
the default values for each parameter will automatically inherit those from that input. This
represents the <i>initial model fit</i>. Proceed to the next tab: <span
style='background-color: lightgrey; color: black; padding: 2px 4px; border-radius: 3px;
font-family: monospace; font-weight: bold;'>Optimize variogram parameters</span>.
</span>

<span style='font-size:18px; line-height:1; margin:0; padding:0; display:block;'>
<b><i>Optimize variogram parameters</i></b></span>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
Variogram model parameter optimization is used to compute the best-fit values for the model
defined within the <span style='background-color: lightgrey; color: black; padding: 2px 4px;
border-radius: 3px; font-family: monospace; font-weight: bold;'>Variogram model</span> tab.
First, click the <span style='background-color: #00bcd4; padding: 2px 4px; border-radius: 3px;
font-family: monospace; font-weight: bold; color: white;'>Update parameter bounds</span> button
to populate the window with the model-specific parameters. Each parameter will have four
options: <b>`Min`</b>, <b>`Value`</b>, <b>`Max`</b>, and <b>`Vary`</b>. This provides
instructions to the built-in optimization algorithm on how to use each variable. When
<b>`Vary`</b> is unchecked, <b>`Value`</b> defines the constant value for that particular
parameter. Otherwise, <b>`Value`</b> is the initial value at the start of the optimization
algorithm. <b>`Min`</b> and <b>`Max`</b> correspond to the minimum and maximum values that the
optimizer is allowed to evaluate, respectively.
</span>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
The <b>`Optimization kwargs`</b> text area represents a formattable dictionary used to inform
the `lmfit`-based optimizer. The optimizer uses the least-squares minimizer, so all `lmfit`
keyword arguments specific to `method="least_squares"` are valid here. This allows for a
dictionary like below to be entered into this textfield that will be passed into the optimizer:
</span>

<pre style='background-color: #f4f4f4; padding: 4px 6px; border-radius: 3px; margin: 2px 0;
font-family: monospace; font-size: 14px; line-height: 1.2;'>
{{
    'max_nfev': 500,
    'ftol': 1e-8,
    'xtol': 1e-8,
    'gtol': 1e-8,
    'diff_step': 1e-10,
    'jac': '3-point'
}}</pre>

<span style='font-size:14px; line-height:1.5; margin:0; padding:0; display:block;'>
After setting up all of the parameter definitions and optimizer keyword arguments, click
<span style='background-color: #337ab7; padding: 2px 4px; border-radius: 3px;
font-family: monospace; font-weight: bold; color: white;'>Optimize variogram parameters</span>
to compute the best-fit variogram parameters. This will produce a <span
style='background-color:red; color:white;'><b>red line</b></span> representing the updated
modeled variogram using the optimized parameters. An additional report will provide a printout
of the best-fit parameters, with bold variables indicating those that were actually minimized
and not held constant. If the GUI instance was stored as a variable, then the parameters can
also be extracted as a dictionary via <span style='color: #4EC9B0; background-color: white;
padding: 2px 4px; border-radius: 3px; font-family: monospace; font-weight: bold;'>
VariogramGUI.optimized_parameters</span>.
</span>

<span style='font-size:18px; line-height:1; margin:0; padding:0; display:block;'>
<b><i>See also</i></b></span>

Documentation for the <span style='color: #4EC9B0; background-color: white; padding: 2px 4px;
border-radius: 3px; font-family: monospace; font-weight: bold;'>
echopop.variogram.Variogram</span> class.
"""

    # Convert to ipywidgets
    INSTRUCTIONS_TAB = ipw.VBox([md_to_HTML(instructions_text)])

    # Return
    return INSTRUCTIONS_TAB


# ==================================================================================================
# VariogramGUI class
# ------------------


class VariogramGUI:

    def __init__(self, data, lag_resolution, n_lags, coordinates, variogram_parameters: dict = {}):

        # Initialize extension
        hv.extension("bokeh")

        # Store data/metadata attributes
        self.data = data

        # Initialize container for `Variogram` object
        self.vgm = None

        # Store the initialized values
        self.initialization_args = dict(
            lag_resolution=lag_resolution, n_lags=n_lags, coordinates=coordinates
        )

        # Store input variogram parameters that will be used as defaults
        self.variogram_parameters = variogram_parameters

        # Initialize the `lmfit.Parameters` and parameter widgets attributes
        self.parameters_lmfit = Parameters()

        # Initialize the input widgets for GUI initialization
        self.initialization_widgets = {}

        # Initialize the parameter widget dictionary for dynamic parameter mapping
        self.parameter_widgets = {}
        self.optimization_parameter_widgets = {}
        self.optimization_parameter_output = None

        # Initialize the plotting area
        self.plot_output = ipw.Output(
            layout=ipw.Layout(border="1px solid black", width="800px", height="500px")
        )
        # ---- Plot state management
        self.plot_layers = {"empirical": None, "theoretical": None, "optimized": None}
        self.empirical_results = {}
        self.theoretical_results = {}
        self.optimization_results = {}

        # Initialize the status pane
        self.status_pane = ipw.HTML(value="", layout=ipw.Layout(width="600px"))

        # Set up the GUI
        self._setup_gui()

    def _setup_gui(self):
        """
        Generate the tabs and overall display for the GUI. These will create both the general
        tabs for each component as well as storing each tab as individual attributes for
        debugging purposes.
        """

        # Instructions tab [TAB 1]
        self.instructions_tab = instructions_tab(self)

        # Initialization tab [TAB 2]
        self.initialization_var_tab = self.initialize_tab()

        # Empirical variogram tab [TAB 3]
        self.empirical_var_tab = self.empirical_tab()

        # Initial variogram model tab [TAB 4]
        self.theoretical_var_tab = self.theoretical_tab()

        # Optimization tab [TAB 5]
        self.optimization_var_tab = self.optimization_tab()

        # Create tabs widget
        control_tabs = ipw.Accordion(
            children=[
                self.instructions_tab,
                self.initialization_var_tab,
                self.empirical_var_tab,
                self.theoretical_var_tab,
                self.optimization_var_tab,
            ]
        )

        # Set which section is open by default
        control_tabs.selected_index = None

        # Update the tab names
        control_tabs.titles = (
            "Instructions",
            "Initialize variogram",
            "Empirical variogram",
            "Variogram model",
            "Optimize variogram parameters",
        )

        # MAIN
        self.tabs = ipw.HBox(
            [
                # ipw.VBox([control_tabs, self.status_pane], layout=ipw.Layout(width="450px")),
                # self.plot_output
                control_tabs,
                ipw.VBox([self.plot_output, self.status_pane]),
            ],
            layout=ipw.Layout(width="100%"),
        )

    @property
    def optimized_parameters(self):
        """
        Return the best-fit, optimized variogram parameters
        """
        return self.optimization_results["optimized_parameters"]

    def _ipython_display_(self):
        """
        Integrate GUI into IPython without needing to create a separate property or 'get'
        method. This enables the usage of just `VariogramGUI(...)` to instantiate a session.
        """
        display(self.tabs)

    def _clear_plot(self):
        """
        Clear the plotting pane
        """
        # Clear the output
        with self.plot_output:
            clear_output(wait=True)

        # Reset all plot layers
        for key in self.plot_layers:
            self.plot_layers[key] = None

    def _clear_downstream_layers(self, from_layer):
        """
        Clear all layers after a given layer in the hierarchy.
        """

        # Get the layer order
        layer_order = list(self.plot_layers.keys())

        # Clear downstream
        if from_layer in layer_order:
            start_idx = layer_order.index(from_layer) + 1
            for layer in layer_order[start_idx:]:
                self.plot_layers[layer] = None

    def _update_plot_display(self):
        """
        Update the plot displayed based on the current plot state.
        """

        with self.plot_output:
            # Clear the plot
            clear_output(wait=True)

            # Get all non-None plot layers in order
            active_plots = [plot for plot in self.plot_layers.values() if plot is not None]

            # Update plot, if any
            if active_plots:
                # ---- Update overlay
                overlay = hv.Overlay(active_plots).opts(
                    hv.opts.Overlay(legend_position="bottom_right")
                )
                # ---- Create dynamic mapper
                self.overlay_dmap = hv.DynamicMap(lambda: overlay)
                # ---- Display
                display(self.overlay_dmap)

    def _update_model_parameters(self, change):
        """
        Dynamically update the displayed variogram model parameter widgets.
        """

        # Get the new model name
        model_name = change["new"]

        # Clear all preexisting parameter widgets
        self.parameter_widgets.clear()

        # Update status panel
        md_to_HTML(
            ipywidget_obj=self.status_pane,
            type="status",
            text=("**Status:** Updating variogram model parameters..."),
        )

        # Refresh the widgets
        with self.parameter_output:
            # ---- Clear
            clear_output(wait=True)
            try:
                # ---- Map the model function
                model = VARIOGRAM_MODEL_PARAMETER_MAP.get(model_name)
                # ---- Get the function signature
                args, _ = get_variogram_arguments(model)
                # ---- Store the model function name
                self.variogram_model = model
                # ---- Extract the parameters
                params = list(p for p in dict(args) if p != "distance_lags")
                # ---- Update the widgets
                for p in params:
                    # ---- Get defaults
                    defaults = DEFAULT_VARIOGRAM_PARAMETERS.get(p, {})
                    # ---- Grab user values
                    user = self.variogram_parameters.get(p, {})
                    # ---- Combine
                    param_info = {**defaults, **user}
                    # ---- Missing parameter? Skip.
                    if not param_info:
                        continue
                    # ---- Generate label
                    param_label = md_to_HTML(
                        f"<p style='line-height: 0.25;'><b>{param_info.get('name', p)}</b></p>"
                    )
                    # ---- Boolean -> Checkbox
                    if isinstance(param_info["value"], bool):
                        param_widget = ipw.Checkbox(
                            value=param_info["value"], layout=ipw.Layout(width="70px")
                        )
                    # ---- Float -> BoundedFloatText
                    else:
                        param_widget = ipw.BoundedFloatText(
                            value=param_info["value"],
                            min=param_info.get("min"),
                            max=param_info.get("max"),
                            step=param_info.get("step"),
                            layout=ipw.Layout(width="70px"),
                        )
                    # ---- Update the parameter widget state
                    self.parameter_widgets[p] = param_widget
                    # ---- Consolidate the full widget
                    widget = ipw.VBox([param_label, param_widget])
                    # ---- Display
                    display(widget)
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="success",
                    text=f"**Status:** Updated {len(params)} parameters for {model_name}.",
                )

                # ---- Port changes
                def create_update_function(param_name):
                    def update_optimization_value(change):
                        if (
                            hasattr(self, "optimization_parameter_widgets")
                            and param_name in self.optimization_parameter_widgets
                        ):
                            self.optimization_parameter_widgets[param_name][
                                "value"
                            ].value = change.new

                    return update_optimization_value

                # ---- Use the factory function to create proper closure
                update_func = create_update_function(p)
                param_widget.observe(update_func, names="value")
            except Exception as e:
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=(
                        f"**Error:** Variogram model could not be computed due to the following "
                        f"error: {str(e)}."
                        f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                    ),
                )

    def _update_optimization_parameters(self):
        """
        Create optimization parameter widgets based on the current theoretical model parameters.
        """

        with self.optimization_parameter_output:
            clear_output(wait=True)

            if not hasattr(self, "parameter_widgets") or not self.parameter_widgets:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=(
                        "**Error:** No theroetical variogram model parameters available. Please "
                        "select a model in the previous tab."
                    ),
                )
                return

            # Update status panel
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="status",
                text=(
                    "**Status:** Updating variogram parameter options available for optimization..."
                ),
            )

            # Iterate through the parameters
            try:
                for param_name, param_widget in self.parameter_widgets.items():
                    # ---- Get the current value from the theoretical tab widget
                    current_value = param_widget.value
                    # ---- Skip boolean parameters for optimization
                    if isinstance(current_value, bool):
                        continue
                    # ---- Get the default parameter info
                    defaults = DEFAULT_VARIOGRAM_PARAMETERS.get(param_name, {})
                    # ---- Create label
                    param_label = md_to_HTML(
                        f"<p style='line-height: 0.25;'><b>{defaults.get('name', param_name)}"
                        f"</b></p>"
                    )
                    # ---- Min value input
                    min_input = ipw.BoundedFloatText(
                        value=defaults.get("min", 0.0),
                        min=-1e10,
                        max=1e10,
                        step=defaults.get("step", 0.1),
                        layout=ipw.Layout(width="75px"),
                    )
                    # ---- Value input (current value from theoretical tab)
                    value_input = ipw.BoundedFloatText(
                        value=current_value if not isinstance(current_value, bool) else 1.0,
                        min=-1e10,
                        max=1e10,
                        step=defaults.get("step", 0.1),
                        layout=ipw.Layout(width="75px"),
                    )
                    # ---- Max value input
                    max_input = ipw.BoundedFloatText(
                        value=defaults.get("max", 100.0),
                        min=-1e10,
                        max=1e10,
                        step=defaults.get("step", 0.1),
                        layout=ipw.Layout(width="75px"),
                    )
                    # ---- Vary checkbox
                    vary_checkbox = ipw.Checkbox(
                        value=defaults.get("vary", True),
                        description="Vary",
                        layout=ipw.Layout(width="70px"),
                        indent=False,
                        disabled=False,
                    )

                    # ---- Function to toggle editability based on vary checkbox
                    def create_toggle_function(min_w, max_w):
                        def toggle_editability(change):
                            min_w.disabled = not change.new
                            max_w.disabled = not change.new

                        return toggle_editability

                    # ---- Set initial state and attach observer
                    toggle_func = create_toggle_function(min_input, max_input)
                    vary_checkbox.observe(toggle_func, names="value")
                    min_input.disabled = not vary_checkbox.value
                    max_input.disabled = not vary_checkbox.value
                    # ---- Store widgets
                    self.optimization_parameter_widgets[param_name] = {
                        "min": min_input,
                        "value": value_input,
                        "max": max_input,
                        "vary": vary_checkbox,
                    }
                    # ---- Create parameter row
                    param_row = ipw.HBox(
                        [
                            ipw.VBox([ipw.HTML(value="<small>Min</small>"), min_input]),
                            ipw.VBox([ipw.HTML(value="<small>Value</small>"), value_input]),
                            ipw.VBox([ipw.HTML(value="<small>Max</small>"), max_input]),
                            ipw.VBox([ipw.HTML(value="<small>Vary</small>"), vary_checkbox]),
                        ],
                        layout=ipw.Layout(width="100%"),
                    )
                    # ---- Complete the parameter section
                    display(ipw.VBox([param_label, param_row]))
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="success",
                    text=(
                        f"**Status**: Optimization parameters (*n*={len(self.parameter_widgets)}) "
                        f"for the variogram model "
                        f"({self.theoretical_widgets['model_selection'].value}) successfully "
                        "loaded."
                    ),
                )
            except Exception as e:
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=(
                        f"**Error:** Variogram model could not be optimized due to the following "
                        f"error: {str(e)}."
                        f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                    ),
                )

    def initialize_tab(self):
        """
        Tab for initializing the empirical variogram calculation. This allows for overriding
        the initial inputs for 'n_lags' and 'lag_resolution'. The column names defined via
        'coordinates' remain static.
        """

        # Input widgets
        # ---- Lag resolution
        lag_resolution_input = ipw.BoundedFloatText(
            value=self.initialization_args["lag_resolution"],
            min=1e-10,
            max=1e99,
            step=self.initialization_args["lag_resolution"] * 1e-1,
            layout=ipw.Layout(width="60px"),
        )
        # -------- Create box
        lag_resolution_box = ipw.VBox(
            [
                md_to_HTML("<p style='line-height: 0.25;'><b>Lag resolution</b></p>"),
                lag_resolution_input,
            ]
        )
        # ---- Number of lags
        n_lags_input = ipw.BoundedIntText(
            value=self.initialization_args["n_lags"],
            min=2,
            max=int(1e99),
            step=1,
            layout=ipw.Layout(width="60px"),
        )
        # -------- Create box
        n_lags_box = ipw.VBox(
            [md_to_HTML("<p style='line-height: 0.25;'><b>Number of lags</b></p>"), n_lags_input]
        )
        # ---- Concatenate horizontally
        input_box = ipw.HBox([lag_resolution_box, HSPACER, n_lags_box])

        # Coordinates text pane
        coord_info = md_to_HTML(f"**Coordinate names:** {self.initialization_args['coordinates']}")

        # Update the widgets container
        self.initialization_widgets.update(
            {"lag_resolution": lag_resolution_input, "n_lags": n_lags_input}
        )

        # Initializion button
        initialize_button = ipw.Button(
            description="Initialize variogram",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Set click interaction
        initialize_button.on_click(self.initialize_variogram)

        return ipw.VBox(
            [
                input_box,
                coord_info,
                VSPACER,
                initialize_button,
            ],
            layout=ipw.Layout(margin="0 0 0 0", padding="0 0 0 0"),
        )

    def initialize_variogram(self, button):
        """
        When the variogram initializatio button is clicked, a `Variogram`-class instance will
        be initialized that enables downstream empirical and theoretical variogram computations.
        """

        try:
            # Get entered inputs
            lag_resolution = self.initialization_widgets["lag_resolution"].value
            n_lags = self.initialization_widgets["n_lags"].value

            # Create new `Variogram`-class instance
            self.vgm = Variogram(
                lag_resolution=lag_resolution,
                n_lags=n_lags,
                coordinate_names=self.initialization_args["coordinates"],
            )

            # Reset the plotting state
            self._clear_plot()

            # Clear large result containers
            self.empirical_results.clear()
            self.theoretical_results.clear()
            self.optimization_results.clear()

            # Update the status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="success",
                text=(
                    f"**Status**: Variogram initialized  \nLag resolution: {lag_resolution}  \n"
                    f"Number of lags: {n_lags}"
                ),
            )
        except Exception as e:
            # Update the status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    f"**Error:** Failed to initialize the variogram due to the following "
                    f"error: {str(e)}. "
                    f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )

    def empirical_tab(self):
        """
        Tab for calculating and visualizing the empirical variogram based on user inputs.
        """

        # Input widgets
        # ---- Variable selection
        variable_select = ipw.Dropdown(
            options=list(VARIABLE_MAP.keys()),
            value="Biomass density (kg nmi⁻²)",
            layout=ipw.Layout(width="250px"),
        )
        # -------- Create box
        variable_select_box = ipw.VBox(
            [md_to_HTML("<p style='line-height: 0.25;'><b>Variable</b></p>"), variable_select]
        )
        # ---- Azimuth filter
        # -------- Set up boolean activation first
        azimuth_filter = ipw.Checkbox(
            value=True, indent=False, layout=ipw.Layout(width="30px", margin="0 5px 0 0")
        )
        # -------- Define threshold widget
        azimuth_threshold = ipw.BoundedFloatText(
            value=180.0, min=0.0, max=180.0, step=1.0, layout=ipw.Layout(width="50px")
        )
        # -------- Degree symbol label
        deg_label = ipw.Label("°", layout=ipw.Layout(width="auto", margin="0 0 0 3px"))
        # ---- Group input + degree label tightly
        azimuth_with_units = ipw.HBox(
            [azimuth_threshold, deg_label], layout=ipw.Layout(align_items="center")
        )

        # -------- Create helper function for toggling the float text accessibility
        def toggle_azimuth_threshold(change):
            azimuth_threshold.disabled = not change.new
            deg_label.layout.display = "none" if azimuth_threshold.disabled else ""

        # -------- Set up observer
        azimuth_filter.observe(toggle_azimuth_threshold, names="value")
        # -------- Hard-code default for filter
        azimuth_threshold.disabled = not azimuth_filter.value
        deg_label.layout.display = "none" if azimuth_threshold.disabled else ""
        # -------- Create azimuth filter box
        azimuth_box = ipw.HBox(
            [azimuth_filter, azimuth_with_units],
            layout=ipw.Layout(justify_content="flex-start", align_items="center"),
        )
        # -------- Add title
        azimuth_input_box = ipw.VBox(
            [
                md_to_HTML(
                    "<p style='line-height: 0.25;'><b>Azimuth angle filter</b> "
                    "(<i>check box to apply filter</i>)</p>"
                ),
                azimuth_box,
            ]
        )

        # Store widgets for callback access
        self.empirical_widgets = {
            "variable": variable_select,
            "azimuth_filter": azimuth_filter,
            "azimuth_threshold": azimuth_threshold,
        }

        # Computation button
        computation_button = ipw.Button(
            description="Compute empirical variogram",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Set click interaction
        computation_button.on_click(self.empirical_variogram)

        return ipw.VBox([variable_select_box, azimuth_input_box, VSPACER, computation_button])

    def empirical_variogram(self, button):
        """
        Calculate the empirical variogram using the updated widget inputs.
        """

        # Disable button
        button.disabled = True

        # Extract the empirical variogram parameters
        variable = VARIABLE_MAP[self.empirical_widgets["variable"].value]
        azimuth_filter = self.empirical_widgets["azimuth_filter"].value
        azimuth_threshold = self.empirical_widgets["azimuth_threshold"].value

        try:
            # Check for initialized `Variogram`-class
            if self.vgm is None:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=(
                        "**Error:** Variogram has not been initialized. Initialize the variogram "
                        "in the 'Initialization' tab before proceeding."
                    ),
                )
                return

            # Check for input dataset
            if self.data is None:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=("**Error:** No dataset provided for variogram computation."),
                )

            # Check if variable exists in data
            if variable not in self.data.columns:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=f"**Error:** Variable '{variable}' not found in dataset.",
                )
                return

            # Compute the empirical variogram
            # ---- Update status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="status",
                text=("**Status:** Computing empirical variogram..."),
            )
            # ---- Compute
            self.vgm.calculate_empirical_variogram(
                data=self.data,
                variable=variable,
                azimuth_filter=azimuth_filter,
                azimuth_angle_threshold=azimuth_threshold,
            )

            # Extract results
            lags = self.vgm.lags
            gamma = self.vgm.gamma
            lag_counts = self.vgm.lag_counts

            # Update the empirical results container
            self.empirical_results = {
                "lags": lags,
                "gamma": gamma,
                "lag_counts": lag_counts,
            }

            # Update status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="success",
                text=(
                    f"**Status:** Empirical variogram computation complete. "
                    f"[{self.empirical_widgets['variable'].value}, *n*={len(lags)} lags]"
                ),
            )

            # Render the plot
            # ---- Store empirical variogram plot
            self.plot_layers["empirical"] = self.plot_empirical_variogram()
            # ---- Clear the current plot state
            self._clear_downstream_layers("empirical")
            # ---- Update the plot display
            self._update_plot_display()

        except Exception as e:
            # Update status panel
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    f"**Error:** Empirical variogram could not be computed due to the following "
                    f"error: {str(e)}."
                    f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )
        finally:
            # Re-enable button
            button.disabled = False

    def plot_empirical_variogram(self):
        """
        Plot the empirical variogram. This is a scatter plot at each defined lag distance where
        the size of each point scales with the lag counts at each point.
        """

        # Extract variables
        lags = self.empirical_results["lags"]
        gamma = self.empirical_results["gamma"]
        lag_counts = self.empirical_results["lag_counts"]

        # Create the required DataFrame for plotting
        plot_data = pd.DataFrame(
            {
                "distance": lags,
                "semivariance": gamma,
                "lag_count": lag_counts,
                "lag_number": range(1, len(lags) + 1),
            }
        )

        # Scale point sizes based on lag counts
        if len(lag_counts) > 1 and lag_counts.max() > lag_counts.min():

            def size_scale(x):
                return (((x - x.min()) / float(x.max() - x.min()) + 1) * 4) ** 2

            plot_data["point_size"] = size_scale(lag_counts)
        else:
            plot_data["point_size"] = 10

        # Create axis-style hook
        def style_axes(plot, element):
            for ax in [plot.state.xaxis[0], plot.state.yaxis[0]]:
                ax.axis_label_text_color = "black"
                ax.axis_label_text_font_size = "16pt"
                ax.axis_label_text_font_style = "normal"
                ax.major_label_text_color = "black"
                ax.major_label_text_font_size = "12pt"

        # Hover tooltip
        hover = HoverTool(
            tooltips=[
                ("Lag #", "@lag_number"),
                ("Distance", "@distance{0.000}"),
                ("Semivariance", "@semivariance{0.000}"),
                ("Count", "@lag_count"),
            ]
        )

        # Create the plot
        variogram_scatter = hv.Scatter(
            plot_data,
            kdims=["distance"],
            vdims=["semivariance", "lag_count", "lag_number", "point_size"],
        ).opts(
            size="point_size",
            color="blue",
            line_color="black",
            line_width=1,
            width=760,
            height=460,
            xlabel="Lag distance (𝒉)",
            ylabel="Semivariance (𝜸)",
            show_grid=True,
            tools=[hover, "pan", "wheel_zoom", "box_zoom", "reset"],
            active_tools=["pan"],
            default_tools=[],
            hooks=[style_axes],
        )

        # Return
        return variogram_scatter

    def theoretical_tab(self):
        """
        Tab for calculating and visualizing the initial theoretical variogram based on user inputs.
        """

        # Initialize the Output container
        self.parameter_output = ipw.Output()

        # Input widgets
        # ---- Model selection
        model_selection = ipw.Dropdown(
            value=None,
            options=list(VARIOGRAM_MODEL_PARAMETER_MAP.keys()),
            layout=ipw.Layout(width="200px"),
        )
        # ---- Attach an observer
        model_selection.observe(self._update_model_parameters, names="value")
        # ---- Create box
        model_selection_box = ipw.VBox(
            [
                md_to_HTML(
                    "<p style='line-height: 0.25;'><b>Variogram model</b> (<i>select one</i>)</p>"
                ),
                model_selection,
            ]
        )

        # Store selection widget
        self.theoretical_widgets = {"model_selection": model_selection}

        # Computation button
        computation_button = ipw.Button(
            description="Compute theoretical variogram",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Set click interaction
        computation_button.on_click(self.theoretical_variogram)

        tab_output = ipw.VBox(
            [model_selection_box, self.parameter_output, VSPACER, computation_button]
        )

        return tab_output

    def theoretical_variogram(self, button):
        """
        Compute the theoretical variogram using the updated widget inputs.
        """

        # Disable button
        button.disabled = True

        # Check for initialized `Variogram`-class
        if self.vgm is None:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    "**Error:** Variogram has not been initialized. Initialize the variogram "
                    "in the 'Initialization' tab before proceeding."
                ),
            )
            return

        # Check for input dataset
        if self.data is None:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=("**Error:** No dataset provided for variogram computation."),
            )

        # Check for empirical result
        if not self.empirical_results:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    "**Error:** You must compute the empirical variogram before running the "
                    "initial theoretical variogram model."
                ),
            )
            return

        # Update status pane
        md_to_HTML(
            ipywidget_obj=self.status_pane,
            type="status",
            text=("**Status:** Computing theoretical variogram..."),
        )

        try:
            # Compute the initial variogram model fit
            # ---- Collect parameter widget inputs
            param_values = {p: w.value for p, w in self.parameter_widgets.items()}
            # ---- Compute the theoretical variogram
            gamma_model = compute_variogram(
                model=self.variogram_model,
                distance_lags=self.vgm.lags,
                variogram_parameters=param_values,
            )

            # Update the initial theoretical model results container
            self.theoretical_results = {
                "lags": self.vgm.lags,
                "gamma": gamma_model,
            }

            # Render the plot
            # ---- Store empirical variogram plot
            self.plot_layers["theoretical"] = self.plot_theoretical_variogram()
            # ---- Clear the current plot state
            self._clear_downstream_layers("theoretical")
            # ---- Update the plot display
            self._update_plot_display()

            # Update status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="success",
                text=(
                    f"**Status:** Theoretical variogram successfully computed "
                    f"({self.theoretical_widgets['model_selection'].value})."
                ),
            )

        except Exception as e:
            # Update status panel
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    f"**Error:** Theoretical variogram model could not be computed due to the  "
                    f"following error: {str(e)}."
                    f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )
        finally:
            # Re-enable button
            button.disabled = False

    def plot_theoretical_variogram(self):
        """
        Plot the theoretical variogram. This is a line curve that fits the modeled semivariance
        predictions to the empirical variogram domain.
        """

        # Extract variables
        lags = self.theoretical_results["lags"]
        gamma = self.theoretical_results["gamma"]

        # Create the required DataFrame for plotting
        plot_data = pd.DataFrame(
            {
                "distance": lags,
                "semivariance": gamma,
            }
        )

        # Create line curve
        variogram_curve = hv.Curve(
            plot_data, kdims=["distance"], vdims=["semivariance"], label="Initial model fit"
        ).opts(color="black", line_width=2.5)

        # Return
        return variogram_curve

    def optimization_tab(self):

        # Initialize the optimization parameter widgets containers
        self.optimization_parameter_output = ipw.Output()
        self.optimization_parameter_widgets = {}

        # Force update optimization parameters when a model is already selected
        if hasattr(self, "parameter_widgets") and self.parameter_widgets:
            self._update_optimization_parameters()

        # Input widgets
        # ---- Optimization dictionary kwargs
        optimization_kwargs_text = ipw.Textarea(
            value="{}",
            placeholder="Enter `lmfit` optimization kwargs as a dictionary...",
            layout=ipw.Layout(width="275px", height="200px"),
        )

        update_params_button = ipw.Button(
            description="Update parameter bounds",
            button_style="info",
            layout=ipw.Layout(width="300px"),
        )
        update_params_button.on_click(lambda b: self._update_optimization_parameters())

        # Optimization button
        optimize_button = ipw.Button(
            description="Optimize variogram parameters",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Attach observer
        optimize_button.on_click(self.optimize_variogram)

        # Store widgets
        self.optimization_widgets = {"kwargs": optimization_kwargs_text}

        # Format title for parameters
        parameter_title = ipw.VBox(
            [
                md_to_HTML(
                    "<span style='font-size:18px; line-height:1; margin:0; padding:0; "
                    "display:block;'>"
                    "<b>Parameter boundaries and variation</b></span>"
                    "<span style='font-size:14px; line-height:1.5; margin:0; padding:0; "
                    "display:block;'><i>Click 'Update parameter bounds' to select which "
                    "parameters to optimize</i></span>"
                )
            ],
            layout=ipw.Layout(margin="0", padding="0"),
        )
        # ---- Create box
        parameter_box = ipw.VBox(
            [parameter_title, update_params_button, self.optimization_parameter_output]
        )

        # Format the UI column
        tab_output = ipw.VBox(
            [
                parameter_box,
                # update_params_button,
                # self.optimization_parameter_output,
                VSPACER,
                ipw.VBox(
                    [
                        md_to_HTML("**Optimization kwargs (lmfit dictionary):**"),
                        optimization_kwargs_text,
                    ]
                ),
                VSPACER,
                optimize_button,
            ]
        )

        return tab_output

    def optimize_variogram(self, button):

        # Disable button
        button.disabled = True

        # Check for initialized `Variogram`-class
        if self.vgm is None:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    "**Error:** Variogram has not been initialized. Initialize the variogram "
                    "in the 'Initialization' tab before proceeding."
                ),
            )
            return

        # Check for input dataset
        if self.data is None:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=("**Error:** No dataset provided for variogram computation."),
            )

        # Check for empirical result
        if not self.empirical_results:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    "**Error:** You must compute the empirical variogram before optimizing "
                    "variogram parameters."
                ),
            )
            return

        try:
            # Update status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="status",
                text="**Status:** Starting parameter optimization...",
            )

            # Parse the varied parameter count (must be at least one)
            vary_count = sum(
                1
                for widgets in self.optimization_parameter_widgets.values()
                if widgets["vary"].value
            )
            # ---- Checksum
            if vary_count == 0:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=(
                        "**Error:** No parameters are set to vary. At least one parameter must be "
                        "varied to perform the variogram model optimization."
                    ),
                )
                return

            # Parse the optimization kwargs entry
            try:
                kwargs_text = self.optimization_widgets["kwargs"].value
                # ---- Convert to a literal dictionary
                optimization_kwargs = ast.literal_eval(kwargs_text)
                # ---- Type check
                if not isinstance(optimization_kwargs, dict):
                    raise ValueError("Optimization kwargs must be a dictionary!")
            except Exception as e:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    type="fail",
                    text=f"**Error:** Invalid optimization kwargs: {str(e)}",
                )
                return

            # Initialize the Parameters object
            lmfit_params = Parameters()
            # ---- Iterate through each parameter to set min/value/max/vary
            for param_name, widgets in self.optimization_parameter_widgets.items():
                min_val = widgets["min"].value
                value_val = widgets["value"].value
                max_val = widgets["max"].value
                vary_val = widgets["vary"].value
                # ---- Add to the object
                lmfit_params.add(
                    name=param_name, value=value_val, min=min_val, max=max_val, vary=vary_val
                )
            # ---- Store
            self.parameters_lmfit = lmfit_params

            # Perform optimization using existing variogram fitting method
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="status",
                text="**Status:** Running optimization...",
            )
            optimization_result = self.vgm.fit_variogram_model(
                model=self.variogram_model,
                model_parameters=lmfit_params,
                optimizer_kwargs=optimization_kwargs,
            )

            # Compute the updated fit
            gamma_optimized = compute_variogram(
                model=self.variogram_model,
                distance_lags=self.vgm.lags,
                variogram_parameters=optimization_result,
            )

            # Store the results
            self.optimization_results = {
                "optimized_parameters": optimization_result,
                "gamma_optimized": gamma_optimized,
            }

            # Store the optimized plot
            self.plot_layers["optimized"] = self.plot_optimized_variogram()
            # ---- Update the plot display
            self._update_plot_display()

            # # Update the status pane
            # ---- Organize the results
            params = self.vgm.variogram_params_optimized
            # ---- Get variation status
            varied_params = [k for k in self.parameters_lmfit if self.parameters_lmfit[k].vary]
            # ---- Format text
            variable_names = [
                (
                    f"**{DEFAULT_VARIOGRAM_PARAMETERS[p]['short_name']}: {v:.5f}**"
                    if p in varied_params
                    else f"{DEFAULT_VARIOGRAM_PARAMETERS[p]['short_name']}: {v:.5f}"
                )
                for p, v in params.items()
            ]
            # ---- Add top-line text
            full_text = [
                "<div style='font-size:18px; font-weight:bold;'>Optimized parameters "
                "<span style='font-weight:normal; font-style:italic;'>"
                "(bold indicates adjusted variables)</span></div>"
            ] + variable_names
            # ---- Final formatting
            text = "  \n".join(full_text)
            # ---- Send to status pane
            md_to_HTML(ipywidget_obj=self.status_pane, text=text)
        except Exception as e:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                type="fail",
                text=(
                    f"**Error:** Optimization failed: {str(e)}  \n"
                    f"**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )
        finally:
            # Re-enable button
            button.disabled = False

    def plot_optimized_variogram(self):
        """
        Plot the optimized theoretical variogram.
        """

        # Create DataFrame
        plot_data = pd.DataFrame(
            {
                "distance": self.vgm.lags,
                "semivariance": self.optimization_results["gamma_optimized"],
            }
        )

        # Create line curve
        variogram_curve = hv.Curve(
            plot_data, kdims=["distance"], vdims=["semivariance"], label="Optimized model fit"
        ).opts(color="red", line_width=2.5)

        return variogram_curve
