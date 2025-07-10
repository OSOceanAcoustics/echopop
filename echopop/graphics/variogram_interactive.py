import warnings
from typing import Optional, Union

import ipywidgets as widgets
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from IPython import get_ipython
from IPython.display import clear_output, display
from matplotlib.patches import Patch

from ..spatial.variogram import (
    empirical_variogram,
    get_variogram_arguments,
    initialize_initial_optimization_values,
    initialize_optimization_config,
    optimize_variogram,
    variogram,
)
from ..utils.validate_dict import (
    VariogramBase,
    VariogramEmpirical,
    VariogramInitial,
    VariogramOptimize,
)

####################################################################################################
# Variogram model dropdown widget map
VARIOGRAM_MODEL_PARAMETER_MAP = {
    "Bessel-exponential": ["bessel", "exponential"],
    "Bessel": "bessel",
    "Exponential": "exponential",
    "Gaussian": "gaussian",
    "Linear": "linear",
    "Sinc": "sinc",
    "Spherical": "spherical",
    "Bessel-Gaussian": ["bessel", "gaussian"],
    "Cosine-exponential": ["cosine", "exponential"],
    "Cosine-Gaussian": ["cosine", "gaussian"],
    "Exponential-linear": ["exponential", "linear"],
    "Gaussian-linear": ["gaussian", "linear"],
}

# Variogram model parameter map
THEORETICAL_VARIOGRAM_PARAMETER_MAP = {
    "nugget": dict(name="Nugget", widget="entry", step=0.01, string_format="0.3f"),
    "sill": dict(name="Sill", widget="entry", step=0.01, string_format="0.3f"),
    "correlation_range": dict(
        name="Correlation range", widget="entry", step=0.0001, string_format="0.4f"
    ),
    "decay_power": dict(name="Decay power", widget="entry", step=0.05, string_format="0.2f"),
    "hole_effect_range": dict(
        name="Hole effect range", widget="entry", step=0.0001, string_format="0.3f"
    ),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox", step=np.nan),
}

# Empirical variogram general settings map
EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {"n_lags": 30, "lag_resolution": 0.002, "azimuth_range": 360.0}

# Optimization
OPTIMIZATION_PARAMETER_MAP = {
    "Maximum function evaluations": dict(
        name="max_fun_evaluations", widget="entry_int", value=500, step=1
    ),
    "Cost function tolerance": dict(
        name="cost_fun_tolerance", widget="entry_float", value=1e-6, step=1e-7
    ),
    "Solution tolerance": dict(
        name="solution_tolerance", widget="entry_float", value=1e-6, step=1e-5
    ),
    "Gradient tolerance": dict(
        name="gradient_tolerance", widget="entry_float", value=1e-4, step=1e-5
    ),
    "Finite differences step size": dict(
        name="finite_step_size", widget="entry_float", value=1e-8, step=1e-9
    ),
    "Trust Region solver method": dict(
        name="trust_region_solver", widget="entry_str", value="exact", step=["exact", "base"]
    ),
    "Characteristic x-scale": dict(
        name="x_scale", widget="entry_hybrid", value="jacobian", value_float=1.0, step=0.1
    ),
    "Jacobian approximation method": dict(
        name="jacobian_approx", widget="entry_str", value="central", step=["central", "forward"]
    ),
}


####################################################################################################
# EMPIRICAL WIDGETS
####################################################################################################
def plot_empirical_variogram(
    lags: np.ndarray,
    gamma_h: np.ndarray,
    lag_counts: np.ndarray,
    patch: Optional[Union[list, Patch]] = None,
):
    """
    Empirical variogram plot panel
    """

    # Generate point-size scaling lambda function
    def size_scale(x):
        return (((x - x.min()) / float(x.max() - x.min()) + 1) * 10) ** 2

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 5))  # Figure and axis initialization
    scatter = ax.scatter(
        lags,
        gamma_h,
        s=size_scale(lag_counts),
        edgecolors="black",
        label="Empirical variogram [γ(h)]",
    )  # scatter plot
    # ---- Patch the legend if needed
    if patch is not None:
        if not isinstance(patch, list):
            ax.legend(handles=[scatter, patch], loc="lower right")
        else:
            ax.legend(handles=[scatter, patch[0], patch[1]], loc="lower right")
    else:
        ax.legend(loc="lower right")
    ax.set_xlabel("Lag distance [h]")  # x-axis label
    ax.set_ylabel("Semivariance [γ]")  # y-axis label
    ax.grid(True)  # create axis grid lines
    ax.set_axisbelow(True)  # set axis grid lines below points
    plt.tight_layout()  # adjust padding around plotting area

    # Generate list of annotation points
    annotation_points = ax.annotate(
        "",
        xy=(0, 0),
        xytext=(-50, -50),
        textcoords="offset points",
        bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
        arrowprops=dict(arrowstyle="->", color="black", lw=1, linestyle="-"),
    )
    # ---- Set values to be invisible
    annotation_points.set_visible(False)

    # Helper function for revealing invisible annotations
    def reveal_annotation(ind):
        position = scatter.get_offsets()[ind["ind"][0]]  # get scatter point positions
        annotation_points.xy = position  # set the annotation coordinates
        # ---- Create text label
        text = f"h={position[0]:.2f}\nγ={position[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}"
        annotation_points.set_text(text)  # assign text to annotation point
        annotation_points.get_bbox_patch().set_alpha(1.0)  # set alpha transparency for text box

    # Helper function for target annotation
    def on_hover(event):
        vis = annotation_points.get_visible()  # detect the visible annotation point
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            # ---- If the annotation intersects with the cursor hover point
            if cont:
                reveal_annotation(ind)
                annotation_points.set_visible(True)
                # ---- Draw
                fig.canvas.draw_idle()
            else:
                if vis:  # points that aren't 'highlighted' by the cursor
                    annotation_points.set_visible(False)
                    # ---- Draw
                    fig.canvas.draw_idle()

    # Create event connection
    fig.canvas.mpl_connect("motion_notify_event", on_hover)

    # Clean up the canvas
    fig.canvas.toolbar_visible = False
    fig.canvas.header_visible = False
    fig.canvas.footer_visible = False
    fig.canvas.resizable = False

    # Return figure and axis
    return fig, ax


####################################################################################################
# EMPIRICAL WIDGETS
####################################################################################################
def empirical_variogram_widget(
    transect_data: pd.DataFrame,
    settings_dict: dict,
    loading_label: widgets.Label,
    plot_output_widget: widgets.Output,
):
    """
    Empirical variogram widgets and interactive elements
    """
    # Pre-allocate style and layout values
    style = {"description_width": "150px"}
    layout = {"width": "400px"}

    # Create general settings value entry widgets
    # ---- Number of lags
    n_lags_entry = widgets.IntText(
        value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
        description="Number of lags",
        step=1,
        style=style,
        layout=layout,
    )
    # ---- Lag (increment) resolution
    lag_resolution_entry = widgets.FloatText(
        value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"],
        description="Lag resolution",
        step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"] / 2,
        style=style,
        layout=layout,
    )
    # ---- Azimuth range
    azimuth_range_entry = widgets.FloatText(
        value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"],
        description="Azimuth Range",
        step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"] / 360,
        style=style,
        layout=layout,
    )

    # Add variable dropdown widget
    empirical_variogram_dropdown = widgets.Dropdown(
        options=list(["Biomass (biomass density)", "Abundance (number density)"]),
        value="Biomass (biomass density)",
        description="Transect variable",
        disabled=False,
        style={"description_width": "150px"},
        layout={"width": "400px"},
    )

    # Add empirical variogram button
    compute_empirical_variogram_button = widgets.Button(
        description="Compute Empirical Variogram",
        layout={"width": "405px"},
        style={"font_weight": "bold"},
    )

    # Initialize `settings_dict` analysis variable
    settings_dict.update({"variable": "biomass_density"})

    # Define updating function based on dropdown selection
    def update_settings_dict(change):
        settings_dict.update(
            {
                "variable": (
                    "biomass_density"
                    if empirical_variogram_dropdown.value == "Biomass (biomass density)"
                    else "number_density"
                ),
            }
        )

    # Initialize the value dictionary
    # Validate the entries and format them accordingly
    # ---- Lag resolution, number of lags
    valid_base_params = VariogramBase.create(
        **{"n_lags": n_lags_entry.value, "lag_resolution": lag_resolution_entry.value}
    )
    # ---- Azimuth range
    valid_empirical_params = VariogramEmpirical.create(
        **{"azimuth_range": azimuth_range_entry.value}
    )

    # Initialize the value dictionary
    general_settings = {
        "n_lags": valid_base_params["n_lags"],
        "lag_resolution": valid_base_params["lag_resolution"],
        "azimuth_range": valid_empirical_params["azimuth_range"],
        "variable": empirical_variogram_dropdown.value,
    }
    # general_settings = {
    #     "n_lags": n_lags_entry.value,
    #     "lag_resolution": lag_resolution_entry.value,
    #     "azimuth_range": azimuth_range_entry.value,
    #     "variable": empirical_variogram_dropdown.value,
    # }

    # Helper function for copmuting the empirical variogram
    def compute_empirical_variogram(
        transect_data: pd.DataFrame, general_settings: dict, settings_dict: dict
    ):
        # ---- Attach the update function to the dropdown value change event
        empirical_variogram_dropdown.observe(update_settings_dict, names="value")
        # ---- Compute the distance lags
        distance_lags = (
            np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
        )
        # ---- Compute the maximum range
        max_range = general_settings["n_lags"] * general_settings["lag_resolution"]
        # ---- Update the general settings for computing the empirical variogram
        general_settings.update(
            {"force_lag_zero": True, "distance_lags": distance_lags, "range": max_range}
        )
        # ---- Compute the empirical variogram
        lags, gamma_h, lag_counts, _ = empirical_variogram(
            transect_data, general_settings, settings_dict
        )
        # ---- Return the outputs
        return lags, gamma_h, lag_counts, general_settings

    # Helper function defining what happens when `Compute empirical variogram` button is pressed
    # ---- Initialize empirical variogram extraction dictionary
    empirical_params = {}

    def click_empirical_variogram_button(button, empirical_params):
        # ---- Show the loading indicator
        loading_label.value = "Computing the empirical variogram, please wait..."
        # ---- Compute the empirical variogram
        lags, gamma_h, lag_counts, _ = compute_empirical_variogram(
            transect_data, general_settings, settings_dict
        )
        # ---- Update `empirical_params`
        empirical_params["lags"] = lags
        empirical_params["gamma_h"] = gamma_h
        empirical_params["lag_counts"] = lag_counts
        # ---- Clear the previous plot
        plot_output_widget.clear_output()
        # ---- Plot the new empirical variogram
        with plot_output_widget:
            fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
            plt.show(fig)
        # ---- Hide the loading indicator
        loading_label.value = ""

    # Attach the button click event to the handler function
    compute_empirical_variogram_button.on_click(
        lambda button: click_empirical_variogram_button(button, empirical_params)
    )

    # Create helper function that updates the current entry values
    def update_general_settings(change):
        # ---- Number of lags
        general_settings["n_lags"] = n_lags_entry.value
        # ---- Lag (increment) resolution
        general_settings["lag_resolution"] = lag_resolution_entry.value
        # ---- Azimuth range
        general_settings["azimuth_range"] = azimuth_range_entry.value

    # Helper function for batch observer registration
    def register_observers(widgets, handler, names="value"):
        for widget in widgets:
            widget.observe(handler, names=names)

    # Register the widget observers
    register_observers(
        [n_lags_entry, lag_resolution_entry, azimuth_range_entry], update_general_settings
    )

    # Vertically concatenate/bundle widgets
    widget_box = widgets.VBox(
        [
            n_lags_entry,
            lag_resolution_entry,
            azimuth_range_entry,
            empirical_variogram_dropdown,
            compute_empirical_variogram_button,
        ]
    )

    # Return the widget and settings dictionary
    return widget_box, general_settings, empirical_params


####################################################################################################
# THEORETICAL WIDGETS
####################################################################################################
def compute_theoretical_variogram(lags: np.ndarray, parameters: dict):
    """
    Compute the theoretical variogram
    """

    # Get the current input values
    arg_dict = {
        key: value.value if isinstance(value, (widgets.FloatText, widgets.Checkbox)) else value
        for key, value in parameters.items()
    }

    # Compute variogram
    return variogram(lags, arg_dict)


def plot_theoretical_variogram(fig, ax, empirical_params: dict, parameters: dict):
    """
    Plot the theoretical variogram
    """
    lags = empirical_params["lags"]
    theoretical_variogram = compute_theoretical_variogram(lags, parameters)

    ax.plot(lags, theoretical_variogram, linewidth=3.0, c="black")
    # plt.plot(lags, theoretical_variogram)
    fig.canvas.draw_idle()


def get_variogram_defaults(variogram_parameters: dict, argument: str):
    """
    Get default variogram parameter values for the particular model
    """
    DEFAULT_PARAMETERS = {
        "nugget": 0.0,
        "sill": 1.0,
        "hole_effect_range": 0.0,
        "decay_power": 1.5,
        "enhance_semivariance": False,
    }
    # ---- get variogram config
    if argument in variogram_parameters.keys():
        return variogram_parameters[argument]
    else:
        return DEFAULT_PARAMETERS[argument]


def theoretical_variogram_widgets(
    empirical_params: dict, variogram_parameters: dict, plot_output_widget: widgets.Output
):

    # Pre-allocate style and layout values
    style = {"description_width": "150px"}
    layout = {"width": "400px"}

    # Dropdown selection with possible model names
    dropdown_variogram_model = widgets.Dropdown(
        options=list(VARIOGRAM_MODEL_PARAMETER_MAP.keys()),
        value="Bessel-exponential",
        description="Variogram model",
        disabled=False,
        style={"description_width": "153px"},
        layout={"width": "403px"},
    )

    # Create theoretical model computation button
    compute_theoretical_variogram_button = widgets.Button(
        description="Compute Theoretical Variogram",
        layout={"width": "405px"},
        style={"font_weight": "bold"},
    )

    # Dynamically populate the model parameter field(s)
    # ---- Initialize the value entry dictionary
    variogram_settings = {}
    # ---- Initialize the Output widget
    variogram_model_widgets = widgets.Output()

    # Helper function for dynamically updating the parameter input widgets
    def dynamic_parameters(change):
        with variogram_model_widgets:
            # ---- Clear previous values
            clear_output(wait=True)
            # ---- Get the dropdown selection
            model_selection = change["new"]
            # ---- Retrieve the model argument used for `variogram(...)`
            model_name = VARIOGRAM_MODEL_PARAMETER_MAP[model_selection]
            # ---- Get argument names for this particular function
            function_arguments, _ = get_variogram_arguments(model_name)
            # ---- Model arguments (convert to dictionary)
            args_dict = dict(function_arguments)
            # ---- Clear `variogram_settings`
            variogram_settings.clear()
            # ---- Clear existing widgets, if any
            variogram_model_widgets.clear_output()
            # ---- Add the model name to the dictionary
            variogram_settings["model"] = model_name
            # ---- Generate widgets for each parameter
            for param in args_dict.values():
                # ---- Get the parameter name
                parameter = param.name
                # ---- Search reference dictionary for widget parameterization
                if parameter != "distance_lags":
                    # ---- Get default parameterization
                    default_value = get_variogram_defaults(variogram_parameters, parameter)
                    # ---- Get the function argument name
                    arg_name = THEORETICAL_VARIOGRAM_PARAMETER_MAP[parameter]["name"]
                    # ---- Get widget-type
                    widget_type = THEORETICAL_VARIOGRAM_PARAMETER_MAP[parameter]["widget"]
                    # ---- Get step size
                    step_size = THEORETICAL_VARIOGRAM_PARAMETER_MAP[parameter]["step"]
                    # ---- Create the appropriate widget
                    if widget_type == "checkbox":
                        checkbox = widgets.Checkbox(
                            value=default_value,
                            description=arg_name,
                            style={"description_width": "initial"},
                        )
                        # ---- Add checkbox entry to `variogram_settings`
                        variogram_settings[parameter] = checkbox
                        # ---- Display the result
                        display(checkbox)
                    elif widget_type == "entry":
                        entry = widgets.FloatText(
                            value=default_value,
                            description=arg_name,
                            style=style,
                            layout=layout,
                            step=step_size,
                        )
                        # ---- Add floattext entry to `variogram_settings`
                        variogram_settings[parameter] = entry
                        # ---- Display the result
                        display(entry)

    # Attach the observer function to the dropdown menu
    dropdown_variogram_model.observe(dynamic_parameters, names="value")

    # Register the updated parameters
    dynamic_parameters({"new": dropdown_variogram_model.value})

    def on_compute_button_clicked(button):
        # parameters = variogram_widgets
        parameters = variogram_settings
        with plot_output_widget:
            # out.clear_output(wait=True)  # Clear previous plot
            clear_output(wait=True)
            patch_initial = Patch(color="black", linewidth=1, label="Initial theoretical variogram")
            fig, ax = plot_empirical_variogram(**empirical_params, patch=patch_initial)
            plot_theoretical_variogram(fig, ax, empirical_params, parameters)
            fig.canvas.toolbar_visible = False
            fig.canvas.header_visible = False
            fig.canvas.footer_visible = False
            fig.canvas.resizable = False
            plt.show(fig)

    # Update plot
    compute_theoretical_variogram_button.on_click(on_compute_button_clicked)

    widget_box = widgets.VBox(
        [dropdown_variogram_model, variogram_model_widgets, compute_theoretical_variogram_button]
    )

    return widget_box, variogram_settings


####################################################################################################
# OPTIMIZE WIDGETS
####################################################################################################
def optimize_variogram_widgets(
    variogram_settings: dict,
    empirical_params: dict,
    general_settings: dict,
    plot_output_widget: widgets.Output,
):
    """
    Variogram optimization widgets
    """

    # Create button 'Optimize variogram parameters' button
    optimize_variogram_button = widgets.Button(
        description="Optimize variogram parameters",
        layout={"width": "405px"},
        style={"font_weight": "bold"},
    )

    # Create HTML widgets to display results
    html_widget = widgets.HTML()

    # Initialize dictionary for tracking argument values
    optimization_args = {}
    float_widgets = {}
    checkbox_widgets = {}

    # Create a container for the widgets
    widget_entries = []

    # Define layout and style
    opt_layout = {"width": "400px"}
    opt_style = {"description_width": "210px"}

    # Helper function to update values dictionary
    def update_values(change):
        widget = change["owner"]
        description = widget.description
        for desc, config in OPTIMIZATION_PARAMETER_MAP.items():
            if desc == description:
                widget_type = config["widget"]
                name = config["name"]

                if widget_type == "entry_int" or widget_type == "entry_float":
                    optimization_args[name] = widget.value
                elif widget_type == "entry_str":
                    optimization_args[name] = widget.value
                elif widget_type == "entry_hybrid":
                    checkbox_value = checkbox_widgets[description].value
                    if checkbox_value:
                        optimization_args[name] = "jacobian"
                    else:
                        optimization_args[name] = float_widgets[description].value

    # Create widgets based on OPTIMIZATION_PARAMETER_MAP
    for description, config in OPTIMIZATION_PARAMETER_MAP.items():
        widget_type = config["widget"]
        name = config["name"]
        initial_value = VariogramOptimize.create(**{name: config["value"]})[name]
        step = config["step"]

        if widget_type == "entry_int":
            widget = widgets.IntText(
                value=initial_value,
                description=description,
                step=step,
                layout=opt_layout,
                style=opt_style,
            )
        elif widget_type == "entry_float":
            widget = widgets.FloatText(
                value=initial_value,
                description=description,
                step=step,
                layout=opt_layout,
                style=opt_style,
            )
        elif widget_type == "entry_str":
            widget = widgets.Dropdown(
                options=step,
                value=initial_value,
                description=description,
                layout=opt_layout,
                style=opt_style,
            )
        elif widget_type == "entry_hybrid":
            float_widget = widgets.FloatText(
                value=config["value_float"],
                description=description,
                layout=opt_layout,
                style=opt_style,
                step=step,
                disabled=True,
            )
            checkbox_widget = widgets.Checkbox(
                description="Use Jacobian for characteristic x-scaling",
                value=True,
                layout=opt_layout,
                style={"description_width": "120px"},
            )

            # Function to toggle between float widget and disabled state
            def toggle_float_widget(change):
                if change.new:
                    optimization_args["x_scale"] = "jacobian"
                    float_widget.disabled = True
                else:
                    float_widget.disabled = False
                    optimization_args["x_scale"] = float(float_widget.value)

            checkbox_widget.observe(toggle_float_widget, "value")

            # Initial state based on default value
            if initial_value == "jacobian":
                float_widget.disabled = True

            widget = widgets.VBox([checkbox_widget, float_widget])

        # Attach observer to each widget to update values dictionary
        widget.observe(update_values, names="value")
        optimization_args[name] = initial_value
        widget_entries.append(widget)

    def compute_optimal_variogram(best_fit_parameters, lags, model_def):

        # Add model name to `best_fit_parameters`
        # ---- Create copy
        best_fit_parameters_copy = best_fit_parameters.copy()
        # ---- Append model
        best_fit_parameters_copy.update({"model": model_def})

        # Compute variogram
        return variogram(lags, best_fit_parameters_copy)

    def plot_optimal_variogram(fig, ax, best_fit_variogram, empirical_params, model):

        # Get lags
        lags = empirical_params["lags"]

        # Compute the optimal variogram
        optimal_variogram = compute_optimal_variogram(best_fit_variogram, lags, model)

        ax.plot(lags, optimal_variogram, linewidth=3.0, c="red", label="Optimal variogram")
        fig.canvas.draw_idle()

    best_fit_variogram_capture = {}

    def on_optimize_variogram_button_clicked(button, variogram_settings, general_settings):
        html_widget.value = ""

        # Convert the widget dictionary into a normal dictionary
        parameters_dict = {
            key: (
                {"value": value.value}
                if isinstance(value, (widgets.FloatText, widgets.Checkbox))
                else value
            )
            for key, value in variogram_settings.items()
        }

        # Extract the model name (this is stored for later)
        model = parameters_dict["model"]
        # ---- Remove from the parameters dictionary
        del parameters_dict["model"]

        # Initialize the starting values
        fit_parameters = initialize_initial_optimization_values(
            VariogramInitial.create(parameters_dict), VariogramBase.create(**general_settings)
        )

        # Get the fit parameters
        # fit_parameters = list(parameters_dict.keys())

        # Create a list of lower bound tuples
        # ---- Repeat zero's
        # lower_bound_values = np.zeros(len(fit_parameters))
        # # ---- Convert to tuples list
        # lower_bounds = [w for w in zip(fit_parameters, lower_bound_values) if None not in w]

        # # Create a tuples list for the upper bounds
        # # ---- Repeat INF
        # upper_bound_values = np.repeat(np.inf, len(fit_parameters))
        # # ---- Convert to tuples list
        # upper_bounds = [w for w in zip(fit_parameters, upper_bound_values) if None not in w]

        # # Create a tuples list for the initial values
        # initial_values = list(parameters_dict.items())

        # Create a refreshed parameter dictionary that includes both the model name and range
        # ---- Create copy
        variogram_parameters = parameters_dict.copy()
        # ---- Update the dictionary
        variogram_parameters.update({"model": model, "range": general_settings["range"]})

        # Create the correctly formatted optimization settings
        # optimization_settings = create_optimization_options(
        #     fit_parameters,
        #     variogram_parameters,
        #     initial_values=initial_values,
        #     lower_bounds=lower_bounds,
        #     upper_bounds=upper_bounds,
        #     model=variogram_parameters["model"],
        #     **optimization_args,
        # )
        optimization_config = initialize_optimization_config(optimization_args)

        # # Validate all variogram-related inputs
        # validate_variogram_parameters(variogram_parameters, fit_parameters, optimization_settings)
        optimization_settings = {"parameters": fit_parameters, "config": optimization_config}

        # Find the best-fit parameters
        best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
            empirical_params["lag_counts"],
            empirical_params["lags"],
            empirical_params["gamma_h"],
            optimization_settings,
            **variogram_parameters,
        )

        # Update plot
        initial_parameters = variogram_settings
        with plot_output_widget:
            clear_output(wait=True)
            patch_initial = Patch(color="black", linewidth=1, label="Initial theoretical variogram")
            patch_optim = Patch(color="red", linewidth=1, label="Optimized variogram")
            fig, ax = plot_empirical_variogram(
                **empirical_params, patch=[patch_initial, patch_optim]
            )
            plot_theoretical_variogram(fig, ax, empirical_params, initial_parameters)
            plot_optimal_variogram(fig, ax, best_fit_variogram, empirical_params, model)
            fig.canvas.toolbar_visible = False
            fig.canvas.header_visible = False
            fig.canvas.footer_visible = False
            fig.canvas.resizable = False
            plt.show(fig)

        # Update the label widgets with the results
        best_fit_variogram_str = "<br>".join(
            f"<b>{THEORETICAL_VARIOGRAM_PARAMETER_MAP[key]['name']}</b>: "
            f"{value:{THEORETICAL_VARIOGRAM_PARAMETER_MAP[key]['string_format']}}"
            for key, value in best_fit_variogram.items()
        )

        # Create an HTML widget
        html_widget.value = (
            f"<div style='font-size: 18px; margin-bottom: -20px;'><b>Best fit variogram:</b></div>"
            f"<br> {best_fit_variogram_str}"
        )

        # Update results values
        best_fit_variogram_capture.update(
            {
                "model_fit": best_fit_variogram,
                "initial_fit": initial_fit[2],
                "optim_fit": optimized_fit[2],
            }
        )
        # best_fit_variogram.update({"initial_fit": initial_fit[2], "optim_fit": optimized_fit[2]})
        # best_fit_variogram["initial_fit"] = initial_fit[2]
        # best_fit_variogram["optim_fit"] = optimized_fit[2]

    optimize_variogram_button.on_click(
        lambda b: on_optimize_variogram_button_clicked(b, variogram_settings, general_settings)
    )

    # Display the widgets
    entry_widgets = widgets.VBox(widget_entries)
    return (
        widgets.VBox([entry_widgets, optimize_variogram_button, html_widget]),
        optimization_args,
        best_fit_variogram_capture,
    )


####################################################################################################
# Group operations
####################################################################################################
def stitch_variogram_accordion(settings_widget, variogram_widget, optimization_widget):
    """
    Stitch together all of the widgets together
    """
    # Initialize the accordion tabs
    variogram_accordion = widgets.Accordion(layout=dict(width="450px"))

    # Populate the accordion children
    variogram_accordion.children = [settings_widget, variogram_widget, optimization_widget]

    # Default first tab
    variogram_accordion.selected_index = 0

    # Set the accordion tab titles
    variogram_accordion.set_title(
        0,
        "Empirical variogram",
    )
    variogram_accordion.set_title(
        1,
        "Theoretical variogram",
    )
    variogram_accordion.set_title(
        2,
        "Optimize variogram parameters",
    )

    # Return the accordion
    return variogram_accordion


####################################################################################################
# Create variogram plotting/fitting GUI
####################################################################################################
def variogram_widgets(
    transect_data: pd.DataFrame,
    variogram_parameters: dict,
    settings_dict: dict,
    analysis_dict: dict,
    results_dict: dict,
):
    """
    GUI for semivariogram fitting
    """

    # Suppress matplotlib- and traitlets-related backend warnings
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    # Enable `matplotlib` widget backend
    # get_ipython().run_line_magic("matplotlib", "widget")
    try:
        # ---- Only works in IPython environments
        get_ipython().run_line_magic("matplotlib", "widget")
    except AttributeError:
        # ---- Fallback for non-IPython environments
        plt.switch_backend("module://ipympl.backend_nbagg")

    # Create Output widget for plotting area
    plot_output_widget = widgets.Output(layout=widgets.Layout(height="550px"))

    # Create Label widget for loading message
    loading_label = widgets.Label(value="")

    # Create empirical variogram widgets
    tab_empirical_variogram, general_settings, empirical_params = empirical_variogram_widget(
        transect_data, settings_dict, loading_label, plot_output_widget
    )

    # Create theoretical variogram widgets
    tab_contents_variogram, variogram_settings = theoretical_variogram_widgets(
        empirical_params, variogram_parameters, plot_output_widget
    )

    # Create variogram optimization widgets
    tab_contents_optimization, optimization_settings, best_fit_results = optimize_variogram_widgets(
        variogram_settings, empirical_params, general_settings, plot_output_widget
    )

    # Create Save Fit Button
    save_fit_button = widgets.Button(description="Save fit", style={"font_weight": "bold"})

    # Stitch together all tabs into an accordion widget
    full_accordion = stitch_variogram_accordion(
        tab_empirical_variogram,
        tab_contents_variogram,
        widgets.VBox([tab_contents_optimization, save_fit_button]),
    )

    # Placeholder for storing results
    result_container = {"best_fit": None, "optimization": None, "variogram": None}

    # Function to return results when Apply Fit is clicked
    def save_fit_clicked(button, result_container, analysis_dict, results_dict, settings_dict):
        # When the button is clicked, store the results
        result_container["best_fit"] = best_fit_results
        result_container["optimization"] = optimization_settings
        result_container["variogram"] = variogram_settings

        # Add optimization settings to the analysis settings
        settings_dict.update({"optimization": result_container["optimization"]})

        # Add "partial" results to analysis attribute
        analysis_dict.update(
            {
                "model": result_container["variogram"]["model"],
                "initial_fit": result_container["best_fit"]["initial_fit"],
                "optimized_fit": result_container["best_fit"]["optim_fit"],
            }
        )

        # Save the results
        results_dict.update(
            {
                "model_fit": result_container["best_fit"]["model_fit"],
                "optimization_settings": result_container["optimization"],
                "model": result_container["variogram"]["model"],
            }
        )

    # Attach click event to the button
    save_fit_button.on_click(
        lambda b: save_fit_clicked(b, result_container, analysis_dict, results_dict, settings_dict)
    )

    # Create widget
    gui_widget = widgets.HBox(
        [widgets.VBox([full_accordion]), widgets.VBox([loading_label, plot_output_widget])]
    )

    # Result container FOR DEBUGGING
    gui_widget.results = result_container
    gui_widget.empirical = empirical_params
    # result_container = {
    #     "empirical_params": empirical_params
    # }

    # Create widget
    return gui_widget  # , results_dict
