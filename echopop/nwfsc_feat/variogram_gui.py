import traceback
from typing import Optional

import holoviews as hv
import ipywidgets as ipw
import markdown
import pandas as pd
from bokeh.models import HoverTool
from IPython.display import clear_output, display
from lmfit import Parameters

from echopop.nwfsc_feat.variogram import Variogram
from echopop.nwfsc_feat.variogram_models import get_variogram_arguments, variogram

# ==================================================================================================
# Constants, mappings, and helper functions API
# ---------------------------------------------
# Vertical spacer
VSPACER = ipw.Box(layout=ipw.Layout(height="10px"))

# Variogram variable mapping
VARIABLE_MAP = {
    "Number density (animals nmi^-2)": "number_density",
    "Biomass density (kg nmi^-2)": "biomass_density",
    "NASC (m^2 nmi^-2)": "nasc",
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
        "min": 1e-10,
        "value": 1.0,
        "max": 99999,
        "vary": False,
        "step": 1.0,
    },
    "decay_power": {
        "name": "Decay power (\u03b1)",
        "min": 1e-10,
        "value": 1.0,
        "max": 2.0,
        "vary": False,
        "step": 0.1,
    },
    "enhance_semivariance": {"name": "Enhance semivariance", "value": True},
    "hole_effect_range": {
        "name": "Hole effect range (a\u2095)",
        "min": 0.0,
        "value": 0.0,
        "max": 99999,
        "vary": False,
        "step": 1.0,
    },
    "sill": {
        "name": "Sill (C)",
        "min": 1e-10,
        "value": 1.0,
        "max": 99999,
        "vary": False,
        "step": 1.0,
    },
    "nugget": {
        "name": "Nugget (C\u2080)",
        "min": 0.0,
        "value": 0.0,
        "max": 99999,
        "vary": False,
        "step": 1.0,
    },
    "smoothness_parameter": {
        "name": "Matérn shape parameter (\u03bd)",
        "min": 0.0,
        "value": 0.5,
        "max": 10.0,
        "vary": False,
        "step": 0.1,
    },
    "shape_parameter": {
        "name": "Scale (\u03b2)",
        "min": 1e-10,
        "value": 1.0,
        "max": 100.0,
        "vary": False,
        "step": 1.0,
    },
    "power_exponent": {
        "name": "Power (\u03c9)",
        "min": 1e-10,
        "value": 1.0,
        "max": 2.0 - 1e-10,
        "vary": False,
        "step": 0.1,
    },
}


def md_to_HTML(
    text: str,
    ipywidget_obj: Optional[ipw.HTML] = None,
):
    """
    Helper function that converts text and regular expressions from a markdown format to a
    `ipywidgets.HTML` text object. When an `ipywidgets.HTML` object is supplied, then the text is
    stored directly into the 'value' attribute of that object.
    """

    if ipywidget_obj:
        ipywidget_obj.value = markdown.markdown(text)
    else:
        return ipw.HTML(value=markdown.markdown(text))


# ==================================================================================================
# Instructions tab
# ----------------
instructions_text = """
# Variogram Analysis Interactive GUI

## Overview
This interactive tool allows you to perform comprehensive variogram analysis with the following
steps:

## Workflow

### Tab 2: Object Initialization
- Set the **lag resolution** and **number of lags** for the variogram analysis
- These parameters control the spatial resolution and extent of the analysis

### Tab 3: Compute Empirical Variogram
- Select the variable to analyze
- Configure azimuth filtering options
- Generate an interactive plot showing:
    - X-axis: Lag distance
    - Y-axis: Semivariance (γ)
    - Point size: Varies with lag count
    - Hover tooltips: Lag #, lag counts, and semivariance values

### Tab 4: Fit Theoretical Model
- Choose from available variogram models
- Parameters dynamically update based on selected model
- Black line shows theoretical fit on empirical plot

### Tab 5: Optimization
- Specify optimization parameters as dictionary
- Set parameter bounds and variation flags
- Red line shows optimized fit
- Display best-fit parameters with appropriate rounding

## Usage Notes
- Make sure to load your data before proceeding
- Each tab builds upon the previous steps
- Results can be saved for further analysis
"""

# Convert to ipywidgets
INSTRUCTIONS_TAB = ipw.VBox([md_to_HTML(instructions_text)])

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
        self.status_pane = ipw.HTML(value="")

        # Set up the GUI
        self._setup_gui()

    def _setup_gui(self):
        """
        Generate the tabs and overall display for the GUI. These will create both the general
        tabs for each component as well as storing each tab as individual attributes for
        debugging purposes.
        """

        # Instructions tab [TAB 1]
        self.instructions_tab = INSTRUCTIONS_TAB

        # Initialization tab [TAB 2]
        self.initialization_var_tab = self.initialize_tab()

        # Empirical variogram tab [TAB 3]
        self.empirical_var_tab = self.empirical_tab()

        # Initial variogram model tab [TAB 4]
        self.theoretical_var_tab = self.theoretical_tab()

        # Create tabs widget
        # self.tabs = ipw.Tab(
        control_tabs = ipw.Accordion(
            children=[
                self.instructions_tab,
                self.initialization_var_tab,
                self.empirical_var_tab,
                self.theoretical_var_tab,
            ]
        )

        # Set which section is open by default
        control_tabs.selected_index = None

        # Update the tab names
        # self.tabs.titles = ("Instructions", "Initialize variogram", "Empirical variogram",
        #                     "Variogram model")
        control_tabs.titles = (
            "Instructions",
            "Initialize variogram",
            "Empirical variogram",
            "Variogram model",
        )

        # MAIN
        self.tabs = ipw.HBox(
            [
                ipw.VBox([control_tabs, self.status_pane], layout=ipw.Layout(width="450px")),
                self.plot_output,
            ],
            layout=ipw.Layout(width="100%"),
        )

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
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                text=f"DEBUG: Found {len(active_plots)} active plots",
            )  # ADD THIS

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
                    param_label = md_to_HTML(f"**{param_info.get('name', p)}**")
                    # ---- Boolean -> Checkbox
                    if isinstance(param_info["value"], bool):
                        param_widget = ipw.Checkbox(
                            value=param_info["value"],
                        )
                    # ---- Float -> BoundedFloatText
                    else:
                        param_widget = ipw.BoundedFloatText(
                            value=param_info["value"],
                            min=param_info.get("min"),
                            max=param_info.get("max"),
                            step=param_info.get("step"),
                        )
                    # ---- Update the parameter widget state
                    self.parameter_widgets[p] = param_widget
                    # ---- Consolidate the full widget
                    widget = ipw.VBox([param_label, param_widget], layout=ipw.Layout(width="200px"))
                    # ---- Display
                    display(widget)
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    text=f"**Status:** Updated {len(params)} parameters for {model_name}.",
                )
            except Exception as e:
                # ---- Update status pane
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
                    text=(
                        f"**Error:** Variogram model could not be computed due to the following "
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
            layout=ipw.Layout(width="300px"),
        )
        # ---- Number of lags
        n_lags_input = ipw.BoundedIntText(
            value=self.initialization_args["n_lags"],
            min=2,
            max=int(1e99),
            step=1,
            layout=ipw.Layout(width="300px"),
        )

        # Coordinates text pane
        coord_info = md_to_HTML(f"**Coordinate names:** {self.initialization_args["coordinates"]}")

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
                ipw.VBox([md_to_HTML("**Lag resolution**"), lag_resolution_input]),
                ipw.VBox([md_to_HTML("**Number of lags**"), n_lags_input]),
                coord_info,
                VSPACER,
                initialize_button,
            ]
        )

        # # Format the UI column
        # ui_column = ipw.VBox([
        #     ipw.VBox([md_to_HTML("**Lag resolution**"), lag_resolution_input]),
        #     ipw.VBox([md_to_HTML("**Number of lags**"), n_lags_input]),
        #     coord_info,
        #     VSPACER,
        #     initialize_button,
        # ], layout=ipw.Layout(width="450px"))

        # # Incorporate the plotting area and status pane
        # tab_output = ipw.VBox([
        #     ipw.HBox([ui_column, self.plot_output], layout=ipw.Layout(width="100%")),
        #     self.status_pane
        # ])

        # # Return
        # return tab_output

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
            self.empirical_results.clear()

            # Update the status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                text=(
                    f"**Status**: Variogram initialized  \nLag resolution: {lag_resolution}  \n"
                    f"Number of lags: {n_lags}"
                ),
            )
        except Exception as e:
            # Update the status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
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
            value="Biomass density (kg nmi^-2)",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Azimuth filter
        # -------- Set up boolean activation first
        azimuth_filter = ipw.Checkbox(value=True)
        # -------- Define threshold widget
        azimuth_threshold = ipw.BoundedFloatText(
            value=180.0, min=0.0, max=180.0, step=1.0, layout=ipw.Layout(width="300px")
        )

        # -------- Create helper function for toggling the float text accessibility
        def toggle_azimuth_threshold(change):
            azimuth_threshold.disabled = not change.new

        # -------- Set up observer
        azimuth_filter.observe(toggle_azimuth_threshold, names="value")
        # -------- Hard-code default for filter
        azimuth_threshold.disabled = not azimuth_filter.value

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

        return ipw.VBox(
            [
                md_to_HTML("**Empirical variogram settings**"),
                ipw.VBox([md_to_HTML("**Variable**"), variable_select]),
                ipw.VBox(
                    [
                        md_to_HTML("**Azimuth filter** (*check box to apply filter*)"),
                        ipw.HBox([azimuth_filter, azimuth_threshold]),
                    ]
                ),
                VSPACER,
                computation_button,
            ]
        )

        # Format the UI column
        # ui_column = ipw.VBox([
        #     md_to_HTML("**Empirical variogram settings**"),
        #     ipw.VBox([md_to_HTML("**Variable**"), variable_select]),
        #     ipw.VBox([
        #         md_to_HTML("**Azimuth filter** (*check box to apply filter*)"),
        #         ipw.HBox([azimuth_filter, azimuth_threshold])
        #     ]),
        #     VSPACER,
        #     computation_button
        # ], layout = ipw.Layout(width="450px"))

        # # Incorporate the plotting area and status pane
        # tab_output = ipw.VBox([
        #     ipw.HBox([ui_column, self.plot_output], layout=ipw.Layout(width="100%")),
        #     self.status_pane
        # ])

        # Return
        # return tab_output

    def empirical_variogram(self, button):
        """
        Calculate the empirical variogram using the updated widget inputs.
        """

        try:
            # Check for initialized `Variogram`-class
            if self.vgm is None:
                md_to_HTML(
                    ipywidget_obj=self.status_pane,
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
                    text=("**Error:** No dataset provided for variogram computation."),
                )

            # Extract the empirical variogram parameters
            variable = VARIABLE_MAP[self.empirical_widgets["variable"].value]
            azimuth_filter = self.empirical_widgets["azimuth_filter"].value
            azimuth_threshold = self.empirical_widgets["azimuth_threshold"].value

            # Compute the empirical variogram
            # ---- Update status pane
            md_to_HTML(
                ipywidget_obj=self.status_pane,
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
                text=(
                    f"**Status:** Empirical variogram computation complete. "
                    f"[{self.empirical_widgets["variable"].value}, *n*={len(lags)} lags]"
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
                text=(
                    f"**Error:** Empirical variogram could not be computed due to the following "
                    f"error: {str(e)}."
                    f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )

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
            size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 3) ** 2
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
            # xlabel="Lag distance [<i>h</i>]",
            xlabel=r"Lag distance [$$h$$]",
            ylabel=r"Semivariance [$$γ$$]",
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
            options=list(VARIOGRAM_MODEL_PARAMETER_MAP.keys()),
            value="Exponential",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Attach an observer
        model_selection.observe(self._update_model_parameters, names="value")
        # self._update_model_parameters({"new": model_selection.value})

        # Computation button
        computation_button = ipw.Button(
            description="Compute theoretical variogram",
            button_style="primary",
            layout=ipw.Layout(width="300px"),
        )
        # ---- Set click interaction
        computation_button.on_click(self.theoretical_variogram)

        tab_output = ipw.VBox(
            [
                md_to_HTML("**Theoretical variogram model settings**"),
                ipw.VBox([md_to_HTML("**Variogram model**"), model_selection]),
                self.parameter_output,
                VSPACER,
                computation_button,
            ]
        )

        # # Initialize the model parameter widgets
        # def delayed_init():
        #     time.sleep(0.1)
        #     self._update_model_parameters({"new": model_selection.value})
        # thread = threading.Thread(target=delayed_init)
        # thread.start()

        # # Format the UI column
        # ui_column = ipw.VBox([
        #     md_to_HTML("**Theoretical variogram model settings**"),
        #     ipw.VBox([md_to_HTML("**Variogram model**"), model_selection]),
        #     self.parameter_output,
        #     VSPACER,
        #     computation_button
        # ], layout = ipw.Layout(width="450px"))

        # # Incorporate the plotting area and status pane
        # tab_output = ipw.VBox([
        #     ipw.HBox([ui_column, self.plot_output], layout=ipw.Layout(width="100%")),
        #     self.status_pane
        # ])

        self._update_model_parameters({"new": model_selection.value})
        # FORCE THE INITIAL PARAMETER UPDATE AFTER A TINY DELAY
        # def force_update():
        #     self._update_model_parameters({"new": "Exponential"})

        # # Schedule the update using IPython's event loop
        # model_selection.observe(lambda x: force_update(), names="value")
        # # Trigger it immediately
        # force_update()

        # # Return
        return tab_output

    def theoretical_variogram(self, button):
        """
        Compute the theoretical variogram using the updated widget inputs.
        """

        # Check for initialized `Variogram`-class
        if self.vgm is None:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
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
                text=("**Error:** No dataset provided for variogram computation."),
            )

        # Check for empirical result
        if not self.empirical_results:
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                text=(
                    "**Error:** You must compute the empirical variogram before running the "
                    "initial theoretical variogram model."
                ),
            )
            return

        try:
            # Compute the initial variogram model fit
            # ---- Collect parameter widget inputs
            param_values = {p: w.value for p, w in self.parameter_widgets.items()}
            # ---- Compute the theoretical variogram
            gamma_model = variogram(
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

        except Exception as e:
            # Update status panel
            md_to_HTML(
                ipywidget_obj=self.status_pane,
                text=(
                    f"**Error:** Theoretical variogram model could not be computed due to the  "
                    f"following error: {str(e)}."
                    f"\n  \n**Traceback:**  \n```  \n{traceback.format_exc()}  \n```"
                ),
            )

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
