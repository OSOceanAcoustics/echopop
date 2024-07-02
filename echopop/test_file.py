from echopop import Survey
import ipywidgets
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display
from echopop.spatial.variogram import empirical_variogram
import matplotlib.pyplot as plt
import numpy as np
import copy
from echopop.spatial.transect import edit_transect_columns
from echopop.spatial.projection import transform_geometry
from echopop.spatial.variogram import get_variogram_arguments
from echopop.analysis import variogram_analysis
from echopop.spatial.mesh import griddify_lag_distances
from echopop.spatial.transect import edit_transect_columns
from echopop.spatial.projection import transform_geometry
from echopop.spatial.variogram import empirical_variogram
import os

OPTIMIZATION_PARAMETER_MAP = {
    "Maximum function evaluations": dict(name="max_fun_evaluations", widget="entry_int", value=500, step=1),
    "Cost function tolerance": dict(name="cost_fun_tolerance", widget="entry_float", value=1e-6, step=1e-7),
    "Solution tolerance": dict(name="solution_tolerance", widget="entry_float", value=1e-6, step=1e-5),
    "Gradient tolerance": dict(name="gradient_tolerance", widget="entry_float", value=1e-4, step=1e-5),
    "Finite differences step size": dict(name="finite_step_size", widget="entry_float", value=1e-8, step=1e-9),
    "Trust Region solver method": dict(name="trust_region_solver", widget="entry_str", value="exact", step=["exact", "base"]),
    "Characteristic x-scale": dict(name="x_scale", widget="entry_hybrid", value="jacobian", step=1.0),
    "Jacobian approximation method": dict(name="jacobian_approx", widget="entry_str", value="central", step=["central", "forward"])
}

def optimize_variogram_widget():

    # Initialize dictionary for tracking argument values
    optimization_args = {}
    float_widgets = {}
    checkbox_widgets = {}

    # Create a container for the widgets
    widget_entries = []

    # Define layout and style
    opt_layout = {"width": "420px"}
    opt_style = {"description_width": "200px"}

    # Function to update values dictionary
    def update_values(change):
        widget = change['owner']
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
                        optimization_args[name] = "default"
                    else:
                        optimization_args[name] = float_widgets[description].value

    # Create widgets based on OPTIMIZATION_PARAMETER_MAP
    for description, config in OPTIMIZATION_PARAMETER_MAP.items():
        widget_type = config["widget"]
        initial_value = config["value"]
        step = config["step"]

        if widget_type == "entry_int":
            widget = widgets.IntText(value=initial_value, description=description, step=step, 
                                     layout=opt_layout, style=opt_style)
        elif widget_type == "entry_float":
            widget = widgets.FloatText(value=initial_value, description=description, step=step, 
                                       layout=opt_layout, style=opt_style)
        elif widget_type == "entry_str":
            widget = widgets.Dropdown(options=step, value=initial_value, description=description, 
                                      layout=opt_layout, style=opt_style)
        elif widget_type == "entry_hybrid":
            float_widget = widgets.FloatText(value=step, description=description, layout=opt_layout, 
                                             style=opt_style)
            checkbox_widget = (
                widgets.Checkbox(description="Use Jacobian for characteristic x-scaling", 
                                 value=True, 
                                 layout={"width": "600px"}, 
                                 style={"description_width": "170px"})
            )

            float_widgets[description] = float_widget
            checkbox_widgets[description] = checkbox_widget
            
            # Function to toggle between float widget and disabled state
            def toggle_float_widget(change):
                if change.new:
                    optimization_args["x_scale"] = "jacobian"
                    float_widget.disabled = True
                else:
                    float_widget.disabled = False
                    optimization_args["x_scale"] = float(float_widget.value)
            
            
            # float_widgets[description] = float_widget
            
            # Initial state based on default value
            # if initial_value == "jacobian":
            #     float_widget.disabled = True
            
            widget = widgets.VBox([checkbox_widget, float_widget])

        # Attach observer to each widget to update values dictionary
        checkbox_widget.observe(toggle_float_widget, 'value')
        widget.observe(update_values, names='value')
        optimization_args[config["name"]] = initial_value 

    # Add the created widget to the list
    widget_entries.append(widget)

    # Return the outputs and optimization parameters
    return widget_entries, 


# Create widgets based on OPTIMIZATION_PARAMETER_MAP
for description, config in OPTIMIZATION_PARAMETER_MAP.items():
    widget_type = config["widget"]
    initial_value = config["value"]
    step = config["step"]

    if widget_type == "entry_int":
        widget = widgets.IntText(value=initial_value, description=description, step=step, layout=opt_layout, style=opt_style)
    elif widget_type == "entry_float":
        widget = widgets.FloatText(value=initial_value, description=description, step=step, layout=opt_layout, style=opt_style)
    elif widget_type == "entry_str":
        widget = widgets.Dropdown(options=step, value=initial_value, description=description, layout=opt_layout, style=opt_style)
    elif widget_type == "entry_hybrid":
        float_widget = widgets.FloatText(value=step, description=description, layout=opt_layout, style=opt_style)
        checkbox_widget = widgets.Checkbox(description="Use Jacobian for characteristic x-scaling", value=True, layout={"width": "600px"}, style={"description_width": "170px"})
        
        # Function to toggle between float widget and disabled state
        def toggle_float_widget(change):
            float_widget.disabled = change.new
        
        checkbox_widget.observe(toggle_float_widget, 'value')
        
        # Initial state based on default value
        if initial_value == "jacobian":
            float_widget.disabled = True
        
        widget = widgets.VBox([checkbox_widget, float_widget])

    # Attach observer to each widget to update values dictionary
    widget.observe(update_values, names='value')
    optimization_args[config["name"]] = initial_value 

    # Add the created widget to the list
    widget_entries.append(widget)

# Display the widgets
widgets.VBox(widget_entries)

# Function to update values dictionary
def update_values(change):
    widget = change['owner']
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
                checkbox_value = widget.children[0].value
                if checkbox_value:
                    optimization_args[name] = "default"
                else:
                    optimization_args[name] = widget.children[1].value

def get_optimization_arguments():
    optimization_args = {}


def empirical_variogram_widget():

    # Create general settings value entry widgets
    # ---- Number of lags
    n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
                                   description="Number of lags",
                                   step=1, style=style, layout=layout)
    # ---- Lag (increment) resolution
    lag_resolution_entry =  widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
                                              description='Lag resolution', 
                                              step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
                                              style=style, layout=layout)
    # ---- Azimuth range 
    azimuth_range_entry = widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
                                            description='Azimuth Range', 
                                            step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
                                            style=style, layout=layout)
    
    # Add variable dropdown widget
    empirical_variogram_dropdown = widgets.Dropdown(
        options=list(["Biomass (biomass density)", "Abundance (number density)"]),
        value="Biomass (biomass density)",
        description="Transect variable",
        disabled=False,
        style = {"description_width": "150px"}, 
        layout = {"width": "400px"}
    )  

    # Add empirical variogram button
    compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
                                                        layout={"width": "405px"},
                                                        style={"font_weight": "bold"})    
    
    settings_dict.update({"variable": "biomass_density"})

    def update_settings_dict(change):
        settings_dict.update({
            "variable": (
                "biomass_density" if empirical_variogram_dropdown.value == "Biomass (biomass density)" 
                else "number_density"
            ),
        })

    # Initialize the value dictionary
    general_settings = {
        "n_lags": n_lags_entry.value,
        "lag_resolution": lag_resolution_entry.value,
        "azimuth_range": azimuth_range_entry.value,
        "variable": empirical_variogram_dropdown.value,
    }

    def compute_empirical_variogram(transect_data, general_settings, settings_dict):
        # Attach the update function to the dropdown value change event
        empirical_variogram_dropdown.observe(update_settings_dict, names='value')
        # Compute the distance lags
        distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
        # ---- Compute the maximum range
        max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

        # Update the general settings for computing the empirical variogram
        general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
                                "range": max_range})
        
        # Compute the empirical variogram
        lags, gamma_h, lag_counts, _ = empirical_variogram(
            transect_data, general_settings, settings_dict
        )

        # Return the outputs
        return lags, gamma_h, lag_counts, general_settings
    
    empirical_params = {}
    # # Button-click
    def click_empirical_variogram_button(button, empirical_params):
        # Show the loading indicator
        loading_label.value = "Computing the empirical variogram, please wait..."
        # Compute the empirical variogram
        
        lags, gamma_h, lag_counts, _ = compute_empirical_variogram(transect_data, general_settings, settings_dict)
        # empirical_params = {}
        empirical_params["lags"] = lags
        empirical_params["gamma_h"] = gamma_h
        empirical_params["lag_counts"] = lag_counts
        # Clear the previous plot
        plot_output_widget.clear_output()

        # Plot the new empirical variogram
        with plot_output_widget:
            fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
            plt.show(fig)

        # Hide the loading indicator
        loading_label.value = ""

        # return empirical_params

    # Attach the button click event to the handler function
    compute_empirical_variogram_button.on_click(lambda button: click_empirical_variogram_button(button, empirical_params))


    # compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

    # Create helper function that updates the current entry values
    def update_general_settings(change):
        # ---- Number of lags
        general_settings["n_lags"] = n_lags_entry.value
        # ---- Lag (increment) resolution
        general_settings["lag_resolution"] = lag_resolution_entry.value
        # ---- Azimuth range 
        general_settings["azimuth_range"] = azimuth_range_entry.value

    def register_observers(widgets, handler, names='value'):
        for widget in widgets:
            widget.observe(handler, names=names)

    # Register the widget observers
    register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
                       update_general_settings)
    
    # Vertically concatenate/bundle widgets
    widget_box = widgets.VBox(
        [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
        compute_empirical_variogram_button]
    )

    # Return the widget and settings dictionary
    return widget_box, general_settings, empirical_params


def stitch_variogram_accordion(settings_widget, variogram_widget, optimization_widget):

    # Initialize the accordion tabs
    variogram_accordion = widgets.Accordion(layout=dict(width="450px"))

    # Populate the accordion children
    variogram_accordion.children = [settings_widget, variogram_widget, optimization_widget]

    # Set the accordion tab titles
    variogram_accordion.set_title(0, "Empirical variogram", )
    variogram_accordion.set_title(1, "Theoretical variogram", )
    variogram_accordion.set_title(2, "Optimize variogram parameters", )

    # Return the accordion
    return variogram_accordion

def get_variogram_defaults(variogram_parameters, argument):
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


def theoretical_variogram_widgets():
    dropdown_variogram_model = widgets.Dropdown(
        options=list(VARIOGRAM_MODEL_PARAMETER_MAP.keys()),
        value="Bessel-exponential",
        description="Variogram model",
        disabled=False,
        style = {"description_width": "153px"}, 
        layout = {"width": "403px"}
    )

    # Create theoretical model computation button
    compute_theoretical_variogram_button = widgets.Button(
        description='Compute Theoretical Variogram',
        layout={"width": "405px"},
        style={'font_weight': 'bold'}
    )

    # Dynamically populate the model parameter field(s)
    # ---- Initialize the value entry dictionary
    variogram_settings = {}
    # ---- Initialize the Output widget
    variogram_model_widgets = widgets.Output()

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
                        checkbox = widgets.Checkbox(value=default_value, description=arg_name, 
                                                    style={"description_width": "initial"})
                        # ---- Add checkbox entry to `variogram_settings`
                        variogram_settings[parameter] = checkbox
                        # ---- Display the result
                        display(checkbox)
                    elif widget_type == "entry": 
                        entry = widgets.FloatText(value=default_value, description=arg_name, 
                                                style=style, layout=layout, step=step_size)
                        # ---- Add floattext entry to `variogram_settings`
                        variogram_settings[parameter] = entry
                        # ---- Display the result
                        display(entry)            

    # Attach the observer function to the dropdown menu
    dropdown_variogram_model.observe(dynamic_parameters, names="value")

    # Register the updated parameters
    dynamic_parameters({"new": dropdown_variogram_model.value})

    widget_box = widgets.VBox([
        dropdown_variogram_model,
        variogram_model_widgets,
        compute_theoretical_variogram_button
    ])

    return widget_box, variogram_settings


    def plot_theoretical_variogram(fig, ax, empirical_params, parameters):
        lags = empirical_params["lags"]
        theoretical_variogram = compute_theoretical_variogram(parameters, lags)
        
        ax.plot(lags, theoretical_variogram, linewidth = 3.0, c = "black")
        # plt.plot(lags, theoretical_variogram)
        fig.canvas.draw_idle()

    def on_compute_button_clicked(button):
        # parameters = variogram_widgets
        parameters = variogram_settings
        with plot_output_widget: 
            # out.clear_output(wait=True)  # Clear previous plot
            clear_output(wait=True)
            fig, ax = plot_empirical_variogram(**empirical_params)
            plot_theoretical_variogram(fig, ax, empirical_params, parameters)
            fig.canvas.toolbar_visible = False
            fig.canvas.header_visible = False
            fig.canvas.footer_visible = False
            fig.canvas.resizable = False
            plt.show(fig)



os.chdir("..")
survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
self = survey
survey.transect_analysis()
self = survey
settings_dict = {"variable": "biomass_density",
                 "stratum_name": "stratum_num",
                 "kriging_parameters": {
                    "longitude_reference": -124.78338,
                    "latitude_offset": 45.0,
                    "longitude_offset": -124.78338,
                 }}
transect_input = copy.deepcopy(self.analysis["transect"])
transect_data = edit_transect_columns(transect_input, settings_dict)
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

VARIOGRAM_MODEL_WIDGET_MAP = {
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
    "Gaussian-linear": ["gaussian", "linear"]
}

# Empirical variogram general settings map
EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {
    "n_lags": 30,
    "lag_resolution": 0.002,
    "azimuth_range": 360.0
}

# Variogram model parameter map
THEORETICAL_VARIOGRAM_PARAMETER_MAP = {
    "nugget": dict(name="Nugget", widget="entry", step=0.01),
    "sill": dict(name="Sill", widget="entry", step=0.01),
    "correlation_range": dict(name="Correlation range", widget="entry", step=0.0001),
    "decay_power": dict(name="Decay power", widget="entry", step=0.05),
    "hole_effect_range": dict(name="Hole effect range", widget="entry", step=0.0001),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox", step=np.nan)
}
# %%
# # DEFINE LAYOUT
import os
os.chdir("..")
import warnings
from echopop.spatial.variogram import variogram
from IPython.display import display, clear_output
# from IPython.display import display, HTML
warnings.filterwarnings("ignore", category=DeprecationWarning)
%matplotlib widget
from echopop import Survey
import ipywidgets
import ipywidgets as widgets
import matplotlib.pyplot as plt
from IPython.display import display
from echopop.spatial.variogram import empirical_variogram
import matplotlib.pyplot as plt
import numpy as np
import copy
from echopop.spatial.transect import edit_transect_columns
from echopop.spatial.projection import transform_geometry
from echopop.spatial.variogram import get_variogram_arguments
from echopop.analysis import variogram_analysis
from echopop.spatial.mesh import griddify_lag_distances
from echopop.spatial.transect import edit_transect_columns
from echopop.spatial.projection import transform_geometry
from echopop.spatial.variogram import empirical_variogram
import os
survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
self = survey
survey.transect_analysis()
self = survey
settings_dict = {"variable": "biomass_density",
                 "stratum_name": "stratum_num",
                 "kriging_parameters": {
                    "longitude_reference": -124.78338,
                    "latitude_offset": 45.0,
                    "longitude_offset": -124.78338,
                 }}
transect_input = copy.deepcopy(self.analysis["transect"])
transect_data = edit_transect_columns(transect_input, settings_dict)
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

VARIOGRAM_MODEL_WIDGET_MAP = {
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
    "Gaussian-linear": ["gaussian", "linear"]
}

# Empirical variogram general settings map
EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {
    "n_lags": 30,
    "lag_resolution": 0.002,
    "azimuth_range": 360.0
}

# Variogram model parameter map
THEORETICAL_VARIOGRAM_PARAMETER_MAP = {
    "nugget": dict(name="Nugget", widget="entry", step=0.01),
    "sill": dict(name="Sill", widget="entry", step=0.01),
    "correlation_range": dict(name="Correlation range", widget="entry", step=0.0001),
    "decay_power": dict(name="Decay power", widget="entry", step=0.05),
    "hole_effect_range": dict(name="Hole effect range", widget="entry", step=0.0001),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox", step=np.nan)
}

# Get variogram parameters
variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
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
    "Gaussian-linear": ["gaussian", "linear"]
}

def stitch_variogram_accordion(settings_widget, variogram_widget):

    # Initialize the accordion tabs
    variogram_accordion = widgets.Accordion(layout=dict(width="450px"))

    # Populate the accordion children
    variogram_accordion.children = [settings_widget, variogram_widget]

    # Set the accordion tab titles
    variogram_accordion.set_title(0, "Empirical variogram", )
    variogram_accordion.set_title(1, "Theoretical variogram", )

    # Return the accordion
    return variogram_accordion

def get_variogram_defaults(variogram_parameters, argument):
    DEFAULT_PARAMETERS = {
        "nugget": 0.0,
        "sill": 1.0,
        "hole_effect_range": 0.0,
        "decay_power": 1.5,
        "enhance_semivariance": False,
    }
    # ---- get variogram config
    for argument in DEFAULT_PARAMETERS.keys():
        if argument in variogram_parameters.keys():
            return variogram_parameters[argument]
        else:
            return DEFAULT_PARAMETERS[argument]

def compute_theoretical_variogram(parameters):
    # Get the current input values
    arg_dict = {
        key: value.value if isinstance(value, (widgets.FloatText, widgets.Checkbox)) else value
        for key, value in parameters.items()
    }

    # pass into variogram
    
    distance_lags = lags
    return variogram(distance_lags, arg_dict)

def plot_theoretical_variogram(fig, ax):
    theoretical_variogram = compute_theoretical_variogram(general_settings)
    
    ax.plot(lags, theoretical_variogram, linewidth = 3.0, c = "black")
    # plt.plot(lags, theoretical_variogram)
    fig.canvas.draw_idle()

def click_theoretical_variogram_button(button, variogram_settings, plot_output_widget):

    # Plot/update the plot widget
    with plot_output_widget:
        # Clear 
        clear_output(wait=True)
        # Figure and axis objects
        fig, ax = plot_empirical_variogram()
        # Add the theoretical variogram line to the plot
        plot_theoretical_variogram(fig, ax, variogram_settings)
        # Edit canvas settings
        fig.canvas.toolbar_visible = False
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        fig.canvas.resizable = False
        # Show the plot
        plt.show(fig)

def on_compute_button_clicked(button):
    # parameters = variogram_widgets
    parameters = variogram_settings
    with plot_output_widget: 
        # out.clear_output(wait=True)  # Clear previous plot
        clear_output(wait=True)
        fig, ax = plot_empirical_variogram()
        plot_theoretical_variogram(fig, ax, parameters)
        fig.canvas.toolbar_visible = False
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        fig.canvas.resizable = False
        plt.show(fig)
        # plt.close(fig)

def theoretical_variogram_widget():

    # Create model dropdown wdiget
    dropdown_variogram_model = widgets.Dropdown(
        options=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()),
        value="Bessel-exponential",
        description="Variogram model",
        disabled=False,
        style = {"description_width": "153px"}, 
        layout = {"width": "403px"}
    )    

    # Create theoretical model computation button
    copmute_theoretical_variogram_button = widgets.Button(
        description='Compute Theoretical Variogram',
        layout={"width": "405px"},
        style={'font_weight': 'bold'}
    )

    # Dynamically populate the model parameter field(s)
    # ---- Initialize the value entry dictionary
    variogram_settings = {}
    # ---- Initialize the Output widget
    variogram_model_widgets = widgets.Output()
    # ---- Helper function for updating the parameter fields based on the model dropdown
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
                        checkbox = widgets.Checkbox(value=default_value, description=arg_name, 
                                                    style={"description_width": "initial"})
                        # ---- Add checkbox entry to `variogram_settings`
                        variogram_settings[parameter] = checkbox
                        # ---- Display the result
                        display(checkbox)
                    elif widget_type == "entry": 
                        entry = widgets.FloatText(value=default_value, description=arg_name, 
                                                style=style, layout=layout, step=step_size)
                        # ---- Add floattext entry to `variogram_settings`
                        variogram_settings[parameter] = entry
                        # ---- Display the result
                        display(entry)            

    # Attach the observer function to the dropdown menu
    dropdown_variogram_model.observe(dynamic_parameters, names="value")

    # Register the updated parameters
    dynamic_parameters({"new": dropdown_variogram_model.value})

    # Update plot
    compute_theoretical_variogram_button.on_click(on_compute_button_clicked)

    # Vertically concatenate/bundle widgets
    widget_box = widgets.VBox([
        dropdown_variogram_model,
        variogram_model_widgets,
        copmute_theoretical_variogram_button
    ])

    # Return the widget and settings dictionary
    return widget_box, variogram_settings

style = {'description_width': '150px'}
layout = {'width': '400px'}
def plot_empirical_variogram(lags: np.ndarray, gamma_h: np.ndarray, lag_counts: np.ndarray):
    # Generate point-size scaling lambda function
    size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 10)**2

    # Create the plot
    fig, ax = plt.subplots(figsize=(8, 5)) # Figure and axis initialization
    scatter = ax.scatter(lags, gamma_h, s=size_scale(lag_counts), edgecolors="black") # scatter plot
    ax.set_xlabel("Lag distance (h)") # x-axis label
    ax.set_ylabel("Semivariance (γ)") # y-axis label
    ax.grid(True) # create axis grid lines
    ax.set_axisbelow(True) # set axis grid lines below points
    plt.tight_layout() # adjust padding around plotting area

    # Generate list of annotation points
    annotation_points = ax.annotate("", xy=(0, 0), xytext=(-50, -50), textcoords="offset points",
                                    bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
                                    arrowprops=dict(arrowstyle="->", color="black", lw=1, 
                                                    linestyle="-"))
    # ---- Set values to be invisible 
    annotation_points.set_visible(False)

    # Helper function for revealing invisible annotations
    def reveal_annotation(ind):
        position = scatter.get_offsets()[ind["ind"][0]] # get scatter point positions
        annotation_points.xy = position # set the annotation coordinates
        # ---- Create text label
        text = f"h={position[0]:.2f}\nγ={position[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}" 
        annotation_points.set_text(text) # assign text to annotation point
        annotation_points.get_bbox_patch().set_alpha(1.0) # set alpha transparency for text box

    # Helper function for target annotation
    def on_hover(event):
        vis = annotation_points.get_visible() # detect the visible annotation point
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            # ---- If the annotation intersects with the cursor hover point
            if cont: 
                reveal_annotation(ind)
                annotation_points.set_visible(True)          
                # ---- Draw
                fig.canvas.draw_idle()
            else:
                if vis: # points that aren't 'highlighted' by the cursor
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

# General settings widget
def empirical_variogram_widget():

    # Create general settings value entry widgets
    # ---- Number of lags
    n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
                                   description="Number of lags",
                                   step=1, style=style, layout=layout)
    # ---- Lag (increment) resolution
    lag_resolution_entry =  widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
                                              description='Lag resolution', 
                                              step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
                                              style=style, layout=layout)
    # ---- Azimuth range 
    azimuth_range_entry = widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
                                            description='Azimuth Range', 
                                            step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
                                            style=style, layout=layout)
    
    # Add variable dropdown widget
    empirical_variogram_dropdown = widgets.Dropdown(
        options=list(["Biomass (biomass density)", "Abundance (number density)"]),
        value="Biomass (biomass density)",
        description="Transect variable",
        disabled=False,
        style = {"description_width": "150px"}, 
        layout = {"width": "400px"}
    )  

    # Add empirical variogram button
    compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
                                                        layout={"width": "405px"},
                                                        style={"font_weight": "bold"})    
    
    # Initialize the value dictionary
    general_settings = {
        "n_lags": n_lags_entry.value,
        "lag_resolution": lag_resolution_entry.value,
        "azimuth_range": azimuth_range_entry.value,
        "variable": empirical_variogram_dropdown.value,
    }

    #     # Define the loading label


    def compute_empirical_variogram(transect_data, general_settings, settings_dict):

        # Compute the distance lags
        distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
        # ---- Compute the maximum range
        max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

        # Update the general settings for computing the empirical variogram
        general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
                                "range": max_range})
        
        # Compute the empirical variogram
        lags, gamma_h, lag_counts, _ = empirical_variogram(
            transect_data, general_settings, settings_dict
        )

        # Return the outputs
        return lags, gamma_h, lag_counts, general_settings

    # # Button-click
    def click_empirical_variogram_button(button):
        # Show the loading indicator
        loading_label.value = "Computing the empirical variogram, please wait..."
        # Compute the empirical variogram
        lags, gamma_h, lag_counts, _ = compute_empirical_variogram(transect_data, general_settings, settings_dict)

        # Clear the previous plot
        plot_output_widget.clear_output()

        # Plot the new empirical variogram
        with plot_output_widget:
            fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
            plt.show(fig)

        # Hide the loading indicator
        loading_label.value = ""

    # Attach the button click event to the handler function
    compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

    # Create helper function that updates the current entry values
    def update_general_settings(change):
        # ---- Number of lags
        general_settings["n_lags"] = n_lags_entry.value
        # ---- Lag (increment) resolution
        general_settings["lag_resolution"] = lag_resolution_entry.value
        # ---- Azimuth range 
        general_settings["azimuth_range"] = azimuth_range_entry.value

    def register_observers(widgets, handler, names='value'):
        for widget in widgets:
            widget.observe(handler, names=names)

    # Register the widget observers
    register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
                       update_general_settings)
    
    # Vertically concatenate/bundle widgets
    widget_box = widgets.VBox(
        [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
        compute_empirical_variogram_button]
    )

    # Return the widget and settings dictionary
    return widget_box, general_settings

loading_label = widgets.Label(value="")
plot_output_widget = widgets.Output(layout=widgets.Layout(height="550px"))
tab_empirical_variogram, general_settings = empirical_variogram_widget()
tab_contents_variogram, variogram_settings = theoretical_variogram_widget()

variogram_accordion = widgets.Accordion(layout=dict(width="550px"))
variogram_accordion.children = [tab_empirical_variogram, tab_contents_variogram]
display(widgets.HBox([variogram_accordion, widgets.VBox([loading_label, plot_output_widget])]))
# display(widgets.HBox([tab_empirical_variogram, widgets.VBox([loading_label, plot_output_widget])]))



# def plot_empirical_variogram(lags: np.ndarray, gamma_h: np.ndarray, lag_counts: np.ndarray):
#     # Generate point-size scaling lambda function
#     size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 10)**2

#     # Create the plot
#     fig, ax = plt.subplots(figsize=(8, 5)) # Figure and axis initialization
#     scatter = ax.scatter(lags, gamma_h, s=size_scale(lag_counts), edgecolors="black") # scatter plot
#     ax.set_xlabel("Lag distance (h)") # x-axis label
#     ax.set_ylabel("Semivariance (γ)") # y-axis label
#     ax.grid(True) # create axis grid lines
#     ax.set_axisbelow(True) # set axis grid lines below points
#     plt.tight_layout() # adjust padding around plotting area

#     # Generate list of annotation points
#     annotation_points = ax.annotate("", xy=(0, 0), xytext=(-50, -50), textcoords="offset points",
#                                     bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
#                                     arrowprops=dict(arrowstyle="->", color="black", lw=1, 
#                                                     linestyle="-"))
#     # ---- Set values to be invisible 
#     annotation_points.set_visible(False)

#     # Helper function for revealing invisible annotations
#     def reveal_annotation(ind):
#         position = scatter.get_offsets()[ind["ind"][0]] # get scatter point positions
#         annotation_points.xy = position # set the annotation coordinates
#         # ---- Create text label
#         text = f"h={position[0]:.2f}\nγ={position[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}" 
#         annotation_points.set_text(text) # assign text to annotation point
#         annotation_points.get_bbox_patch().set_alpha(1.0) # set alpha transparency for text box

#     # Helper function for target annotation
#     def on_hover(event):
#         vis = annotation_points.get_visible() # detect the visible annotation point
#         if event.inaxes == ax:
#             cont, ind = scatter.contains(event)
#             # ---- If the annotation intersects with the cursor hover point
#             if cont: 
#                 reveal_annotation(ind)
#                 annotation_points.set_visible(True)          
#                 # ---- Draw
#                 fig.canvas.draw_idle()
#             else:
#                 if vis: # points that aren't 'highlighted' by the cursor
#                     annotation_points.set_visible(False)
#                     # ---- Draw
#                     fig.canvas.draw_idle()

#     # Create event connection
#     fig.canvas.mpl_connect("motion_notify_event", on_hover)

#     # Clean up the canvas
#     fig.canvas.toolbar_visible = False
#     fig.canvas.header_visible = False
#     fig.canvas.footer_visible = False
#     fig.canvas.resizable = False

#     # Return figure and axis
#     return fig, ax

# def compute_empirical_variogram(transect_data, general_settings, settings_dict):

#     # Compute the distance lags
#     distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
#     # ---- Compute the maximum range
#     max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

#     # Update the general settings for computing the empirical variogram
#     general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
#                             "range": max_range})
    
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = empirical_variogram(
#         transect_data, general_settings, settings_dict
#     )

#     # Return the outputs
#     return lags, gamma_h, lag_counts, general_settings

# # Button-click
# def click_empirical_variogram_button(button):
#     # Show the loading indicator
#     loading_label.value = "Computing the empirical variogram, please wait..."
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = compute_empirical_variogram(transect_data, general_settings, settings_dict)

#     # Clear the previous plot
#     plot_output_widget.clear_output()

#     # Plot the new empirical variogram
#     with plot_output_widget:
#         fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
#         plt.show(fig)

#     # Hide the loading indicator
#     loading_label.value = ""

# # # Register the widget observers
# # register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
# #                     update_general_settings)

# # General settings widget
# def empirical_variogram_widget(settings_dict, plot_output_widget, style, layout):

#     # Create general settings value entry widgets
#     # ---- Number of lags
#     n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
#                                    description="Number of lags",
#                                    step=1, style=style, layout=layout)
#     # ---- Lag (increment) resolution
#     lag_resolution_entry =  widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
#                                               description='Lag resolution', 
#                                               step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
#                                               style=style, layout=layout)
#     # ---- Azimuth range 
#     azimuth_range_entry = widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
#                                             description='Azimuth Range', 
#                                             step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
#                                             style=style, layout=layout)
    
   
#     # Add variable dropdown widget
#     empirical_variogram_dropdown = widgets.Dropdown(
#         options=list(["Biomass (biomass density)", "Abundance (number density)"]),
#         value="Biomass (biomass density)",
#         description="Transect variable",
#         disabled=False,
#         style = {"description_width": "150px"}, 
#         layout = {"width": "400px"}
#     )  

#     # Add empirical variogram button
#     compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
#                                                         layout={"width": "405px"},
#                                                         style={"font_weight": "bold"})    
    
#     # Define the loading label
#     loading_label = widgets.Label(value="")
    
#     # Initialize the value dictionary
#     general_settings = {
#         "n_lags": n_lags_entry.value,
#         "lag_resolution": lag_resolution_entry.value,
#         "azimuth_range": azimuth_range_entry.value,
#         "variable": empirical_variogram_dropdown.value,
#     }

#     # Attach the button click event to the handler function
#     compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

#     # Create helper function that updates the current entry values
#     def update_general_settings(change):
#         # ---- Number of lags
#         general_settings["n_lags"] = n_lags_entry.value
#         # ---- Lag (increment) resolution
#         general_settings["lag_resolution"] = lag_resolution_entry.value
#         # ---- Azimuth range 
#         general_settings["azimuth_range"] = azimuth_range_entry.value

#     # Register the widget observers
#     register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
#                        update_general_settings)
    
#     # Vertically concatenate/bundle widgets
#     widget_box = widgets.VBox(
#         [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
#         compute_empirical_variogram_button]
#     )

#     # Return the widget and settings dictionary
#     return widget_box, general_settings

# %matplotlib widget

# import warnings
# from IPython.display import display, HTML
# warnings.filterwarnings("ignore", category=DeprecationWarning)

# %matplotlib widget
# # DEFINE LAYOUT
# style = {'description_width': '150px'}
# layout = {'width': '400px'}

# # Define the loading label
# loading_label = widgets.Label(value="")
# plot_output_widget = widgets.Output(layout=widgets.Layout(height="550px"))

# # Display the button and the plot output widget
# # Empirical variogram widgets
# tab_contents_settings, general_settings = empirical_variogram_widget(settings_dict, plot_output_widget, style, layout)
# # display(compute_empirical_variogram_button, plot_output_widget)

# # Initialize the accordion tabs
# variogram_accordion = widgets.Accordion(layout=dict(width="450px"))
# # Populate the accordion children
# variogram_accordion.children = [tab_contents_settings]
# # Display the tabs
# # tab_contents_settings
# display(widgets.HBox([variogram_accordion, widgets.VBox([loading_label, plot_output_widget])]))

# # Get the stratum name
# stratum_name = self.analysis["settings"]["transect"]["stratum_name"]

# # Get standardization config for kriging 
# standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
# # ---- Get isobath data
# isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# # Get variogram parameters
# variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
# style = {'description_width': '150px'}
# layout = {'width': '400px'}
# widget_settings = dict(style=style, layout=layout)

# # Variogram model dropdown widget map
# VARIOGRAM_MODEL_WIDGET_MAP = {
#     "Bessel-exponential": ["bessel", "exponential"],
#     "Bessel": "bessel",
#     "Exponential": "exponential",
#     "Gaussian": "gaussian",
#     "Linear": "linear",
#     "Sinc": "sinc",
#     "Spherical": "spherical",
#     "Bessel-Gaussian": ["bessel", "gaussian"],
#     "Cosine-exponential": ["cosine", "exponential"],
#     "Cosine-Gaussian": ["cosine", "gaussian"],
#     "Exponential-linear": ["exponential", "linear"],
#     "Gaussian-linear": ["gaussian", "linear"]
# }

# # Empirical variogram general settings map
# EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {
#     "n_lags": 30,
#     "lag_resolution": 0.002,
#     "azimuth_range": 360.0
# }

# # Variogram model parameter map
# THEORETICAL_VARIOGRAM_PARAMETER_MAP = {
#     "nugget": dict(name="Nugget", widget="entry", step=0.01),
#     "sill": dict(name="Sill", widget="entry", step=0.01),
#     "correlation_range": dict(name="Correlation range", widget="entry", step=0.0001),
#     "decay_power": dict(name="Decay power", widget="entry", step=0.05),
#     "hole_effect_range": dict(name="Hole effect range", widget="entry", step=0.0001),
#     "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox", step=np.nan)
# }

# def plot_empirical_variogram(lags: np.ndarray, gamma_h: np.ndarray, lag_counts: np.ndarray):
#     # Generate point-size scaling lambda function
#     size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 10)**2

#     # Create the plot
#     fig, ax = plt.subplots(figsize=(8, 5)) # Figure and axis initialization
#     scatter = ax.scatter(lags, gamma_h, s=size_scale(lag_counts), edgecolors="black") # scatter plot
#     ax.set_xlabel("Lag distance (h)") # x-axis label
#     ax.set_ylabel("Semivariance (γ)") # y-axis label
#     ax.grid(True) # create axis grid lines
#     ax.set_axisbelow(True) # set axis grid lines below points
#     plt.tight_layout() # adjust padding around plotting area

#     # Generate list of annotation points
#     annotation_points = ax.annotate("", xy=(0, 0), xytext=(-50, -50), textcoords="offset points",
#                                     bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
#                                     arrowprops=dict(arrowstyle="->", color="black", lw=1, 
#                                                     linestyle="-"))
#     # ---- Set values to be invisible 
#     annotation_points.set_visible(False)

#     # Helper function for revealing invisible annotations
#     def reveal_annotation(ind):
#         position = scatter.get_offsets()[ind["ind"][0]] # get scatter point positions
#         annotation_points.xy = position # set the annotation coordinates
#         # ---- Create text label
#         text = f"h={position[0]:.2f}\nγ={position[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}" 
#         annotation_points.set_text(text) # assign text to annotation point
#         annotation_points.get_bbox_patch().set_alpha(1.0) # set alpha transparency for text box

#     # Helper function for target annotation
#     def on_hover(event):
#         vis = annotation_points.get_visible() # detect the visible annotation point
#         if event.inaxes == ax:
#             cont, ind = scatter.contains(event)
#             # ---- If the annotation intersects with the cursor hover point
#             if cont: 
#                 reveal_annotation(ind)
#                 annotation_points.set_visible(True)          
#                 # ---- Draw
#                 fig.canvas.draw_idle()
#             else:
#                 if vis: # points that aren't 'highlighted' by the cursor
#                     annotation_points.set_visible(False)
#                     # ---- Draw
#                     fig.canvas.draw_idle()

#     # Create event connection
#     fig.canvas.mpl_connect("motion_notify_event", on_hover)

#     # Clean up the canvas
#     fig.canvas.toolbar_visible = False
#     fig.canvas.header_visible = False
#     fig.canvas.footer_visible = False
#     fig.canvas.resizable = False

#     # Return figure and axis
#     return fig, ax

# def compute_empirical_variogram(transect_data, general_settings, settings_dict):

#     # Compute the distance lags
#     distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
#     # ---- Compute the maximum range
#     max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

#     # Update the general settings for computing the empirical variogram
#     general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
#                             "range": max_range})
    
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = empirical_variogram(
#         transect_data, general_settings, settings_dict
#     )

#     # Return the outputs
#     return lags, gamma_h, lag_counts, general_settings

# # Button-click
# def click_empirical_variogram_button(button):
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = compute_empirical_variogram(transect_data, general_settings, settings_dict)

#     # Clear the previous plot
#     plot_output_widget.clear_output()

#     # Plot the new empirical variogram
#     with plot_output_widget:
#         fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
#         plt.show(fig)

# # Register the widget observers
# register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
#                     update_general_settings)

# # General settings widget
# def empirical_variogram_widget():

#     # Create general settings value entry widgets
#     # ---- Number of lags
#     n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
#                                    description="Number of lags",
#                                    step=1, style=style, layout=layout)
#     # ---- Lag (increment) resolution
#     lag_resolution_entry =  widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
#                                               description='Lag resolution', 
#                                               step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
#                                               style=style, layout=layout)
#     # ---- Azimuth range 
#     azimuth_range_entry = widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
#                                             description='Azimuth Range', 
#                                             step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
#                                             style=style, layout=layout)
    
#     # Add variable dropdown widget
#     empirical_variogram_dropdown = widgets.Dropdown(
#         options=list(["Biomass (biomass density)", "Abundance (number density)"]),
#         value="Biomass (biomass density)",
#         description="Transect variable",
#         disabled=False,
#         style = {"description_width": "150px"}, 
#         layout = {"width": "400px"}
#     )  

#     # Add empirical variogram button
#     compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
#                                                         layout={"width": "405px"},
#                                                         style={"font_weight": "bold"})    
    
#     # Initialize the value dictionary
#     general_settings = {
#         "n_lags": n_lags_entry.value,
#         "lag_resolution": lag_resolution_entry.value,
#         "azimuth_range": azimuth_range_entry.value,
#         "variable": empirical_variogram_dropdown.value,
#     }

#     # Attach the button click event to the handler function
#     compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

#     # Create helper function that updates the current entry values
#     def update_general_settings(change):
#         # ---- Number of lags
#         general_settings["n_lags"] = n_lags_entry.value
#         # ---- Lag (increment) resolution
#         general_settings["lag_resolution"] = lag_resolution_entry.value
#         # ---- Azimuth range 
#         general_settings["azimuth_range"] = azimuth_range_entry.value

#     def register_observers(widgets, handler, names='value'):
#         for widget in widgets:
#             widget.observe(handler, names=names)

#     # Register the widget observers
#     register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
#                        update_general_settings)
    
#     # Vertically concatenate/bundle widgets
#     widget_box = widgets.VBox(
#         [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
#         compute_empirical_variogram_button]
#     )

#     # Return the widget and settings dictionary
#     return widget_box, general_settings

# # General settings widget
# def empirical_variogram_widget(style, layout):

#     # Create general settings value entry widgets
#     # ---- Number of lags
#     n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
#                                    description="Number of lags",
#                                    step=1, style=style, layout=layout)
#     # ---- Lag (increment) resolution
#     lag_resolution_entry =  widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
#                                               description='Lag resolution', 
#                                               step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
#                                               style=style, layout=layout)
#     # ---- Azimuth range 
#     azimuth_range_entry = widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
#                                             description='Azimuth Range', 
#                                             step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
#                                             style=style, layout=layout)
    
#     # Add variable dropdown widget
#     empirical_variogram_dropdown = widgets.Dropdown(
#         options=list(["Biomass (biomass density)", "Abundance (number density)"]),
#         value="Biomass (biomass density)",
#         description="Transect variable",
#         disabled=False,
#         style = {"description_width": "150px"}, 
#         layout = {"width": "400px"}
#     )  

#     # Add empirical variogram button
#     compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
#                                                         layout={"width": "405px"},
#                                                         style={"font_weight": "bold"})    
    
#     # Initialize the value dictionary
#     general_settings = {
#         "n_lags": n_lags_entry.value,
#         "lag_resolution": lag_resolution_entry.value,
#         "azimuth_range": azimuth_range_entry.value,
#     }
    
#     settings_dict.update({
#         "variable": (
#             "biomass" if empirical_variogram_dropdown.value == "Biomass (biomass density)" 
#             else "abundance"
#         ),
#     })

#     # Attach the button click event to the handler function
#     compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

#     # Create helper function that updates the current entry values
#     def update_general_settings(change):
#         # ---- Number of lags
#         general_settings["n_lags"] = n_lags_entry.value
#         # ---- Lag (increment) resolution
#         general_settings["lag_resolution"] = lag_resolution_entry.value
#         # ---- Azimuth range 
#         general_settings["azimuth_range"] = azimuth_range_entry.value

    
#     # Vertically concatenate/bundle widgets
#     widget_box = widgets.VBox(
#         [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
#         compute_empirical_variogram_button]
#     )

#     # Return the widget and settings dictionary
#     return widget_box, general_settings
# ########################
# # Empirical variogram computation function
# def compute_empirical_variogram(transect_data, general_settings, settings_dict):

#     # Compute the distance lags
#     distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
#     # ---- Compute the maximum range
#     max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

#     # Update the general settings for computing the empirical variogram
#     general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, 
#                             "range": max_range})
    
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = empirical_variogram(
#         transect_data, general_settings, settings_dict
#     )

#     # Return the outputs
#     return lags, gamma_h, lag_counts, general_settings

# # Button-click
# def click_empirical_variogram_button(button):
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = compute_empirical_variogram()

#     # Clear the previous plot
#     plot_output_widget.clear_output()

#     # Plot the new empirical variogram
#     with plot_output_widget:
#         fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
#         plt.show(fig)

# ########################
# def empirical_variogram_widget(general_settings, settings_dict, plot_output_widget, style, layout):

#     # Initialize variogram settings for text widgets
#     # ---- Number of lags
#     n_lags_entry = widgets.IntText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["n_lags"],
#                                    description="Number of lags",
#                                    step=1, style=style, layout=layout)
#     # ---- Lag (increment) resolution
#     lag_resolution_entry = (
#         widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"], 
#                           description='Lag resolution', 
#                           step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["lag_resolution"]/2, 
#                           style=style, layout=layout)
#     )
#     # ---- Azimuth range 
#     azimuth_range_entry = (
#         widgets.FloatText(value=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"], 
#                           description='Azimuth Range', 
#                           step=EMPIRICAL_VARIOGRAM_PARAMETER_MAP["azimuth_range"]/360, 
#                           style=style, layout=layout)
#     )
    
#     # Add variable dropdown widget
#     empirical_variogram_dropdown = widgets.Dropdown(
#         options=list(["Biomass (biomass density)", "Abundance (number density)"]),
#         value="Biomass (biomass density)",
#         description="Transect variable",
#         disabled=False,
#         style = {"description_width": "150px"}, 
#         layout = {"width": "400px"}
#     )  

#     # Add empirical variogram button
#     compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
#                                                         layout={"width": "405px"},
#                                                         style={"font_weight": "bold"})    
    


# def compute_empirical_variogram():
#     transect_input = copy.deepcopy(self.analysis["transect"])
#     # ---- Edit the transect data
#     transect_data = edit_transect_columns(transect_input, settings_dict)
#     isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
#     transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)
    
#     # Compute the distance lags
#     distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
#     # Compute the maximum range
#     max_range = general_settings["n_lags"] * general_settings["lag_resolution"]

#     # Update the general settings for computing the empirical variogram
#     general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags, "range": max_range})

#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts, _ = empirical_variogram(
#         transect_data, general_settings, settings_dict
#     )
#     return lags, gamma_h, lag_counts

# def plot_empirical_variogram(lags: np.ndarray, gamma_h: np.ndarray, lag_counts: np.ndarray):
#     # Generate point-size scaling lambda function
#     size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 10)**2

#     # Create the plot
#     fig, ax = plt.subplots(figsize=(8, 5)) # Figure and axis initialization
#     scatter = ax.scatter(lags, gamma_h, s=size_scale(lag_counts), edgecolors="black") # scatter plot
#     ax.set_xlabel("Lag distance (h)") # x-axis label
#     ax.set_ylabel("Semivariance (γ)") # y-axis label
#     ax.grid(True) # create axis grid lines
#     ax.set_axisbelow(True) # set axis grid lines below points
#     plt.tight_layout() # adjust padding around plotting area
#     # Return figure and axis
#     return fig, ax

# def on_compute_button_click(button):
#     # Show the loading indicator
#     loading_label.value = "Computing the empirical variogram, please wait..."
#     # Compute the empirical variogram
#     lags, gamma_h, lag_counts = compute_empirical_variogram()

#     # Clear the previous plot
#     plot_output_widget.clear_output()

#     # Plot the new empirical variogram
#     with plot_output_widget:
#         fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
#         plt.show(fig)

#     # Hide the loading indicator
#     loading_label.value = ""

# # Attach the button click event to the handler function
# compute_button.on_click(on_compute_button_click)

# # Display the button and the plot output widget
# display(compute_button, empirical_variogram_dropdown, loading_label, plot_output_widget)

# def variogram_gui(self):

#     # Get the stratum name
#     stratum_name = self.analysis["settings"]["transect"]["stratum_name"]

#     # Get standardization config for kriging 
#     standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
#     # ---- Get isobath data
#     isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

#     # Get transect data
#     transect_input = copy.deepcopy(self.analysis["transect"])

#     # Generate settings dictionary
#     settings_dict = {
#         "stratum_name": stratum_name,
#         "verbose": False,
#         "kriging_parameters": {
#             "longitude_reference": standardization_parameters["longitude_reference"],
#             "longitude_offset": standardization_parameters["longitude_offset"],
#             "latitude_offset": standardization_parameters["latitude_offset"],
#         },
#     }


#     # Initialize the value dictionary
#     general_settings = {
#         "n_lags": n_lags_entry.value,
#         "lag_resolution": lag_resolution_entry.value,
#         "azimuth_range": azimuth_range_entry.value,
#         "variable": dropdown_variable.value,
#     }

#     # Create helper function that updates the current entry values
#     def update_general_settings(change):
#         # ---- Number of lags
#         general_settings["n_lags"] = n_lags_entry.value
#         # ---- Lag (increment) resolution
#         general_settings["lag_resolution"] = lag_resolution_entry.value
#         # ---- Azimuth range 
#         general_settings["azimuth_range"] = azimuth_range_entry.value

#     # Register the widget observers
#     register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
#                        update_general_settings)
    
#     # React to the "compute empirical variogram button" being clicked
#     lags, gamma_h, lag_counts, general_settings = (
#         compute_empirical_variogram_button.on_click(empirical_variogram_button)
#     )



#     # Generate settings dictionary
#     settings_dict = {
#         "stratum_name": stratum_name,
#         "variable": "biomass" if general_settings["variable"] == "Biomass (biomass density)" else "abundance",
#         "verbose": False,
#         "kriging_parameters": {
#             "longitude_reference": standardization_parameters["longitude_reference"],
#             "longitude_offset": standardization_parameters["longitude_offset"],
#             "latitude_offset": standardization_parameters["latitude_offset"],
#         },
#     }

#     # Create variogram settings dictionary necessayr for preparing transect data
# general_settings = {'n_lags': 30,
#  'lag_resolution': 0.002,
#  'azimuth_range': 360.0,
#  'variable': 'Biomass (biomass density)'}

# # --- Get transect params
# stratum_name = self.analysis["settings"]["transect"]["stratum_name"]
# # --- Get standardization config for kriging 
# standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
# settings_dict = {
#     "stratum_name": stratum_name,
#     "variable": "biomass" if general_settings["variable"] == "Biomass (biomass density)" else "abundance",
#     "verbose": False,
#     "kriging_parameters": {
#         "longitude_reference": standardization_parameters["longitude_reference"],
#         "longitude_offset": standardization_parameters["longitude_offset"],
#         "latitude_offset": standardization_parameters["latitude_offset"],
#     },
# }
# # ---- Create a transect data copy 
# transect_input = copy.deepcopy(self.analysis["transect"])
# # ---- Edit the transect data
# transect_data = edit_transect_columns(transect_input, settings_dict)

# transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

# # Compute distance lags
# distance_lags = np.arange(1, general_settings["n_lags"]) * general_settings["lag_resolution"]
# general_settings.update({"force_lag_zero": True, "distance_lags": distance_lags})
# lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
#     transect_data, general_settings, settings_dict
# )


# # Define the widget map
# VARIOGRAM_MODEL_WIDGET_MAP = {
#     "Bessel-exponential": ["bessel", "exponential"],
#     "Bessel": "bessel",
#     "Exponential": "exponential",
#     "Gaussian": "gaussian",
#     "Linear": "linear",
#     "Sinc": "sinc",
#     "Spherical": "spherical",
#     "Bessel-Gaussian": ["bessel", "gaussian"],
#     "Cosine-exponential": ["cosine", "exponential"],
#     "Cosine-Gaussian": ["cosine", "gaussian"],
#     "Exponential-linear": ["exponential", "linear"],
#     "Gaussian-linear": ["gaussian", "linear"]
# }

# DEFAULT_SETTINGS_PARAMETERS = {
#     "n_lags": 30,
#     "lag_resolution": 0.002,
#     "azimuth_range": 360.0
# }

# def stitch_variogram_accordion(settings_widget, variogram_widget, optimization_widget):

#     # Initialize the accordion tabs
#     variogram_accordion = widgets.Accordion(layout=dict(width="450px"))

#     # Populate the accordion children
#     variogram_accordion.children = [settings_widget, variogram_widget, optimization_widget]

#     # Set the accordion tab titles
#     variogram_accordion.set_title(0, "Empirical variogram", )
#     variogram_accordion.set_title(1, "Theoretical variogram", )
#     variogram_accordion.set_title(2, "Optimize variogram parameters", )

#     # Return the accordion
#     return variogram_accordion

# def plot_empirical_variogram():
#     # Generate point-size scaling lambda function
#     size_scale = lambda x: (((x - x.min()) / float(x.max() - x.min()) + 1) * 10)**2

#     # Create the plot
#     fig, ax = plt.subplots(figsize=(8, 5)) # Figure and axis initialization
#     scatter = ax.scatter(lags, gamma_h, s=size_scale(lag_counts), edgecolors="black") # scatter plot
#     ax.set_xlabel("Lag distance (h)") # x-axis label
#     ax.set_ylabel("Semivariance (γ)") # y-axis label
#     ax.grid(True) # create axis grid lines
#     ax.set_axisbelow(True) # set axis grid lines below points
#     plt.tight_layout() # adjust padding around plotting area

#     # Generate list of annotation points
#     annotation_points = ax.annotate("", xy=(0, 0), xytext=(-50, -50), textcoords="offset points",
#                                     bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
#                                     arrowprops=dict(arrowstyle="->", color="black", lw=1, 
#                                                     linestyle="-"))
#     # ---- Set values to be invisible 
#     annotation_points.set_visible(False)

#     # Helper function for revealing invisible annotations
#     def reveal_annotation(ind):
#         position = scatter.get_offsets()[ind["ind"][0]] # get scatter point positions
#         annotation_points.xy = position # set the annotation coordinates
#         # ---- Create text label
#         text = f"h={position[0]:.2f}\nγ={position[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}" 
#         annotation_points.set_text(text) # assign text to annotation point
#         annotation_points.bbox_patch().set_alpha(1.0) # set alpha transparency for text box

#     # Helper function for target annotation
#     def on_hover(event):
#         vis = annotation_points.get_visible() # detect the visible annotation point
#         if event.inaxes == ax:
#             cont, ind = scatter.contains(event)
#             # ---- If the annotation intersects with the cursor hover point
#             if cont: 
#                 reveal_annotation(ind)
#                 annotation_points.set_visible(True)
#                 fig.canvas.draw_idle()
#             else:
#                 if vis: # points that aren't 'highlighted' by the cursor
#                     annotation_points.set_visible(False)
#                     fig.canvas.draw_idle()

#     # Create event connection
#     fig.canvas.mpl_connect("motion_notify_event", on_hover)

#     # Return figure and axis
#     return fig, ax

# def variogram_plot_widget():

#     # Initialize the Output widget
#     plot_output_widget = widgets.Output(layout=widgets.Layout(height="550px"))

#     # Attach the plot to the Output widget
#     with plot_output_widget:
#         fig, ax = plot_empirical_variogram()
#         plt.show(fig)

# transect_dict =  self.analysis["transect"]
# settings_dict =  self.analysis["settings"]["variogram"]
# isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# # Extract specific variogram parameters
# # ---- Number of lags
# n_lags = variogram_parameters["n_lags"]
# # ---- Lag resolution
# lag_resolution = variogram_parameters["lag_resolution"]

# # Compute the lag distances
# distance_lags = np.arange(1, n_lags) * lag_resolution
# # ---- Add to the `variogram_parameters` dictionary
# variogram_parameters["distance_lags"] = distance_lags
# # ---- Update the max range parameter, if necessary
# max_range = lag_resolution * n_lags


# optimization_parameters["solution_tolerance"] = 1e-4
# optimization_parameters["cost_fun_tolerance"] = 1e-6
# optimization_parameters["gradient_tolerance"] = 1e-4
# optimization_parameters["jacobian_approx"] = "central"
# # Generate the optimization settings dictionary
# optimization_settings = create_optimization_options(
#     fit_parameters,
#     variogram_parameters,
#     initial_values=initial_values,
#     lower_bounds=lower_bounds,
#     upper_bounds=upper_bounds,
#     model=variogram_parameters["model"],
#     **optimization_parameters,
# )

# # Validate all variogram-related inputs
# validate_variogram_parameters(variogram_parameters, fit_parameters, optimization_settings)

# best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
#     lag_counts, lags, gamma_h, variogram_parameters, optimization_settings
# )
# best_fit_variogram

# survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
#                 survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
# survey.transect_analysis()
# survey.fit_variogram(jacobian_approx="central")
# self = survey
# settings_dict = {"variable": "biomass_density",
#                  "stratum_name": "stratum_num",
#                  "kriging_parameters": {
#                     "longitude_reference": -124.78338,
#                     "latitude_offset": 45.0,
#                     "longitude_offset": -124.78338,
#                  }}

# # Prepare the transect data
# # ---- Create a copy of the transect dictionary
# transect_input = copy.deepcopy(self.analysis["transect"])
# # ---- Edit the transect data
# transect_data = edit_transect_columns(transect_input, settings_dict)
# isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
# transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

# # Create a copy of the existing variogram settings
# variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
# n_lags = 30
# lag_resolution = variogram_parameters["lag_resolution"]
# # ---- Number of lags
# variogram_parameters["n_lags"] = 30
# # ---- Azimuth range
# variogram_parameters["azimuth_range"] = 360
# # ---- Force lag-0
# variogram_parameters["force_lag_zero"] = True
# # Compute the lag distances
# distance_lags = np.arange(1, n_lags) * lag_resolution
# # ---- Add to the `variogram_parameters` dictionary
# variogram_parameters["distance_lags"] = distance_lags
# # ---- Update the max range parameter, if necessary
# max_range = lag_resolution * n_lags
# lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
#     transect_data, variogram_parameters, {"variable": "biomass_density"}
# )
# fig1, ax1 = plt.subplots()
# ax1.imshow(np.random.random((5,5)))
# plt.show()
# def plot_variogram():
#     s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
#     x = lags; y = gamma_h

#     fig, ax = plt.subplots(figsize=(8, 6))
#     scatter = ax.scatter(x, y, s=s(lag_counts))
#     ax.set_xlabel('Lag Distance')
#     ax.set_ylabel('Semivariance')
#     ax.set_title('Empirical Variogram')
#     ax.grid(True)
#     plt.tight_layout()
#     return fig

# fig = plot_variogram()

# plot_variogram()
# plt.show()

# def plot_variogram():
#     s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
#     plt.figure(figsize=(6,4))
#     plt.scatter(lags, gamma_h, s=s(lag_counts))
#     plt.xlabel('Lag Distance')
#     plt.ylabel('Semivariance')
#     plt.title('Empirical Variogram')
#     plt.grid(True)
#     plt.tight_layout()

#     # Embed plot into tkinter GUI
#     canvas = FigureCanvasTkAgg(plt.gcf(), master=plot_frame)
#     canvas.draw()
#     canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

#     # Optionally add toolbar for plot navigation
#     toolbar = NavigationToolbar2Tk(canvas, plot_frame)
#     toolbar.update()
#     canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# VARIOGRAM_MODEL_WIDGET_MAP = {
#     "Bessel-exponential": ["bessel", "exponential"],
#     "Bessel": "bessel",
#     "Exponential": "exponential",
#     "Gaussian": "gaussian",
#     "Linear": "linear",
#     "Sinc": "sinc",
#     "Spherical": "spherical",
#     "Bessel-Gaussian": ["bessel", "gaussian"],
#     "Cosine-exponential": ["cosine", "exponential"],
#     "Cosine-Gaussian": ["cosine", "gaussian"],
#     "Exponential-linear": ["exponential", "linear"],
#     "Gaussian-linear": ["gaussian", "linear"]
# }

# optimization_fields = [
#     ("Maximum function evaluations", "max_fun_evaluations"),
#     ("Cost function tolerance", "cost_fun_tolerance"),
#     ("Solution tolerance", "solution_tolerance"),
#     ("Gradient tolerance", "gradient_tolerance"),
#     ("Finite differences step size", "finite_step_size"),
#     ("Trust Region solver method", "trust_region_solver"),
#     ("Characteristic x-scale", "x_scale"),
#     ("Jacobian approximation method", "jacobian_approx")
# ]


# param_key = {
#     "nugget": dict(name="Nugget", widget="entry"),
#     "sill": dict(name="Sill", widget="entry"),
#     "correlation_range": dict(name="Correlation range", widget="entry"),
#     "decay_power": dict(name="Decay power", widget="entry"),
#     "hole_effect_range": dict(name="Hole effect range", widget="entry"),
#     "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox")
# }

# optimization_key = {
#     "max_fun_evaluations": dict(name="Maximum function evaluations", widget="entry"),
#     "cost_fun_tolerance": dict(name="Cost function tolerance", widget="entry"),
#     "solution_tolerance": dict(name="Solution tolerance", widget="entry"),
#     "gradient_tolerance": dict(name="Gradient tolerance", widget="entry"),
#     "finite_step_size": dict(name="Finite differences step size", widget="entry"),
#     "trust_region_solver": dict(name="Trust Region solver method", widget="dropdown"),
#     "x_scale": dict(name="Characteristic x-scale", widget=["checkbox", "entry"]),
#     "jacobian_approx": dict(name="Jacobian approximation method", widget="dropdown"),
# }

# optimization_parameters = {
#     "max_fun_evaluations": max_fun_evaluations,
#     "cost_fun_tolerance": cost_fun_tolerance,
#     "solution_tolerance": solution_tolerance,
#     "gradient_tolerance": gradient_tolerance,
#     "finite_step_size": finite_step_size,
#     "trust_region_solver": trust_region_solver,
#     "x_scale": x_scale,
#     "jacobian_approx": jacobian_approx,
# }

# def get_optimization_defaults(optimization_parameters):
#     # ----


# def on_selection(event):
#     """
#     """
#     selected_model = combobox.get()
#     get_variogram_model(selected_model)

# def get_variogram_model(selected_model = "Bessel-exponential"):
#     """
#     """

#     # Get the corresponding value from the `WIDGET_MAP`
#     model = VARIOGRAM_MODEL_WIDGET_MAP.get(selected_model)

#     # Begin processing
#     if isinstance(model, list):
#         print(f"Processing composite model: {selected_model}")
#     else:
#         print(f"Processing model: {selected_model}")

# # Function to create the optimization parameters frame
# def create_optimization_frame():
#     for field, var_name in optimization_fields:
#         label = ttk.Label(optimization_frame, text=field)
#         label.pack(padx=10, pady=2)
#         if var_name in ["force_lag_zero", "verbose"]:
#             var = tk.BooleanVar(value=True)
#             checkbutton = ttk.Checkbutton(optimization_frame, variable=var)
#             checkbutton.pack(padx=10, pady=2)
#         else:
#             entry = ttk.Entry(optimization_frame)
#             entry.pack(padx=10, pady=2)
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# # Create the main application window
# root = tk.Tk()
# root.title("Fit Variogram GUI")

# # Create a main frame to hold all components
# main_frame = ttk.Frame(root)
# main_frame.pack(fill=tk.BOTH, expand=True)

# # Create a frame for the left side (tabs)
# left_frame = ttk.Frame(main_frame)
# left_frame.grid(row=0, column=0, sticky=tk.NSEW, padx=10, pady=10)

# # Create a frame for the right side (plot panel)
# right_frame = ttk.Frame(main_frame)
# right_frame.grid(row=0, column=1, sticky=tk.NSEW, padx=10, pady=10)

# # Create tabs using ttk.Notebook in the left frame
# notebook = ttk.Notebook(left_frame)
# notebook.pack(fill=tk.BOTH, expand=True)

# # Create frames for each tab
# variogram_frame = ttk.Frame(notebook)
# optimization_frame = ttk.Frame(notebook)

# # Add tabs to the notebook
# notebook.add(variogram_frame, text='Variogram Parameters')
# notebook.add(optimization_frame, text='Optimization Parameters')

# # Dropdown widget
# combobox = ttk.Combobox(left_frame, values=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()), state="readonly", width=20)
# combobox.set("Bessel-exponential")
# combobox.pack(padx=10, pady=5)





# # Function to handle x-scale selection
# def toggle_x_scale():
#     if x_scale_var.get():
#         x_scale_entry.config(state=tk.DISABLED)
#         x_scale_entry.delete(0, tk.END)  # Clear the entry field when disabled
#     else:
#         x_scale_entry.config(state=tk.NORMAL)
        


# # Function to update fit parameters based on selected model
# def update_fit_parameters(event):

#     param_key = {
#         "nugget": "Nugget",
#         "sill": "Sill",
#         "correlation_range": "Correlation range",
#         "decay_power": "Decay power",
#         "hole_effect_range": "Hole effect range"
#     }
#     selected_model = combobox.get()
#     print(f"Selected Model: {selected_model}")

#     # Clear existing widgets in the variogram frame
#     for widget in variogram_frame.winfo_children():
#         widget.destroy()

#     # Generate checkboxes and entry fields based on selected model
#     if selected_model in VARIOGRAM_MODEL_WIDGET_MAP:
#         fit_params = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
#         function_arguments, _ = get_variogram_arguments(fit_params)

#         for param in function_arguments.values():
#             if param.name != "distance_lags":
#                 arg = param.name
#                 # ---- get the name from the key
#                 arg_name = param_key[arg]
#                 var = tk.BooleanVar(value=True)
#                 checkbox = ttk.Checkbutton(variogram_frame, text=arg_name, variable=var)
#                 checkbox.pack(padx=5, pady=5)

#                 # Add an entry field next to each checkbox
#                 entry = ttk.Entry(variogram_frame, width=10)
#                 entry.pack(padx=5, pady=5)
#     else:
#         print(f"No fit parameters defined for model: {selected_model}")

# # Bind selection event to update fit parameters
# combobox.bind("<<ComboboxSelected>>", update_fit_parameters)

# # Function to create the optimization parameters frame
# def create_optimization_frame():
#     for field, var_name in optimization_fields:
#         label = ttk.Label(optimization_frame, text=field)
#         label.pack(padx=10, pady=2)

#         if var_name == "trust_region_solver":
#             combobox = ttk.Combobox(optimization_frame, values=["Exact", "Base"], state="readonly")
#             combobox.pack(padx=10, pady=2)
#         elif var_name == "jacobian_approx":
#             combobox = ttk.Combobox(optimization_frame, values=["Forward", "Central"], state="readonly")
#             combobox.pack(padx=10, pady=2)
#         else:
#             entry = ttk.Entry(optimization_frame)
#             entry.pack(padx=10, pady=2)

# # Initialize optimization frame
# # create_optimization_frame()
# # Create widgets for each optimization parameter
# for field, var_name in optimization_fields:
#     label = ttk.Label(optimization_frame, text=field)
#     label.pack(padx=10, pady=2)

#     if var_name == "x_scale":
#         x_scale_var = tk.BooleanVar(value=True)  # Default checked
#         x_scale_checkbox = ttk.Checkbutton(optimization_frame, text="Use Jacobian", variable=x_scale_var,
#                                            command=toggle_x_scale)
#         x_scale_checkbox.pack(padx=10, pady=2)

#         x_scale_entry = ttk.Entry(optimization_frame, state=tk.DISABLED)
#         x_scale_entry.pack(padx=10, pady=2)

#         # Initial setup based on default value
#         toggle_x_scale()  # Call toggle function initially

#     else:
#         entry = ttk.Entry(optimization_frame)
#         entry.pack(padx=10, pady=2)
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# def plot_variogram():
#     s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
#     x = lags; y = gamma_h
#     # plt.figure(figsize=(8,6))
#     # plt.scatter(x, y, s=s(lag_counts))
#     # plt.xlabel('Lag Distance')
#     # plt.ylabel('Semivariance')
#     # plt.title('Empirical Variogram')
#     # plt.grid(True)
#     # plt.tight_layout()

#     fig, ax = plt.subplots(figsize=(8, 6))
#     scatter = ax.scatter(x, y, s=s(lag_counts))
#     ax.set_xlabel('Lag Distance')
#     ax.set_ylabel('Semivariance')
#     ax.set_title('Empirical Variogram')
#     ax.grid(True)
#     plt.tight_layout()

#     # Embed plot into tkinter GUI
#     canvas = FigureCanvasTkAgg(plt.gcf(), master=right_frame)
#     canvas.draw()
#     canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

#     # Optionally add toolbar for plot navigation
#     toolbar = NavigationToolbar2Tk(canvas, right_frame)
#     toolbar.update()
#     canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

#      # Create a list of annotated points
#     annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
#                             bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
#                             arrowprops=dict(arrowstyle="->", color="black", lw=1, linestyle="-"))
#     annot.set_visible(False)

#     def update_annot(ind):
#         pos = scatter.get_offsets()[ind["ind"][0]]
#         annot.xy = pos
#         text = f"Lag={pos[0]:.2f}\nγ={pos[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}"
#         annot.set_text(text)
#         annot.get_bbox_patch().set_alpha(1.0)

#     def on_hover(event):
#         vis = annot.get_visible()
#         if event.inaxes == ax:
#             cont, ind = scatter.contains(event)
#             if cont:
#                 update_annot(ind)
#                 annot.set_visible(True)
#                 fig.canvas.draw_idle()
#             else:
#                 if vis:
#                     annot.set_visible(False)
#                     fig.canvas.draw_idle()

#     fig.canvas.mpl_connect("motion_notify_event", on_hover)

# # Button to compute empirical variogram
# button_compute = ttk.Button(left_frame, text="Compute Empirical Variogram", command=plot_variogram)
# button_compute.pack(padx=10, pady=10)

# # Start the main event loop
# root.mainloop()