import ipywidgets as widgets
import matplotlib.pyplot as plt
from typing import Dict
from IPython.display import display, clear_output

import pandas as pd
import numpy as np

from echopop.spatial.variogram import (
    get_variogram_arguments, empirical_variogram
)

from echopop.spatial.projection import transform_geometry

from echopop.spatial.transect import edit_transect_columns


# Get the stratum name
stratum_name = self.analysis["settings"]["transect"]["stratum_name"]

# Get standardization config for kriging 
standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
# ---- Get isobath data
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# Get variogram parameters
variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()

# Get transect data
transect_input = copy.deepcopy(self.analysis["transect"])

# Generate settings dictionary

from echopop import Survey

survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
self = survey
survey.transect_analysis()

# Get the stratum name
stratum_name = self.analysis["settings"]["transect"]["stratum_name"]

# Get standardization config for kriging 
standardization_parameters = self.input["statistics"]["kriging"]["model_config"]
# ---- Get isobath data
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# Get variogram parameters
variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()

# Get transect data
transect_input = copy.deepcopy(self.analysis["transect"])


settings_dict = {
    "stratum_name": stratum_name,
    "verbose": False,
    "kriging_parameters": {
        "longitude_reference": standardization_parameters["longitude_reference"],
        "longitude_offset": standardization_parameters["longitude_offset"],
        "latitude_offset": standardization_parameters["latitude_offset"],
    },
}

# Empirical variogram general settings map
EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {
    "n_lags": 30,
    "lag_resolution": 0.002,
    "azimuth_range": 360.0
}

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
EMPIRICAL_VARIOGRAM_PARAMETER_MAP = {
    "n_lags": 30,
    "lag_resolution": 0.002,
    "azimuth_range": 360.0
}

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

# Variogram model parameter map
THEORETICAL_VARIOGRAM_PARAMETER_MAP = {
    "nugget": dict(name="Nugget", widget="entry", step=0.01),
    "sill": dict(name="Sill", widget="entry", step=0.01),
    "correlation_range": dict(name="Correlation range", widget="entry", step=0.0001),
    "decay_power": dict(name="Decay power", widget="entry", step=0.05),
    "hole_effect_range": dict(name="Hole effect range", widget="entry", step=0.0001),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox", step=np.nan)
}

import ipywidgets as widgets
import matplotlib.pyplot as plt
from typing import Dict
from IPython.display import display, clear_output

import pandas as pd
import numpy as np

from echopop.spatial.variogram import (
    get_variogram_arguments, empirical_variogram
)

from echopop.spatial.projection import transform_geometry

from echopop.spatial.transect import edit_transect_columns
style = {'description_width': '150px'}
layout = {'width': '400px'}
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


# Define your dropdown widget
dropdown_variogram_model = widgets.Dropdown(
    options=list(VARIOGRAM_MODEL_PARAMETER_MAP.keys()),
    value="Bessel-exponential",
    description="Variogram model",
    disabled=False,
    style = {"description_width": "153px"}, 
    layout = {"width": "403px"}
)

# Create text entry widgets (initially empty)
text_widgets = {
    param_name: widgets.FloatText(description=param_details['name'], step=param_details['step'])
    if param_details['widget'] == 'entry' else
    widgets.Checkbox(description=param_details['name'])
    for param_name, param_details in THEORETICAL_VARIOGRAM_PARAMETER_MAP.items()
}

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
])

widget_box


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

# Button-click
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

    return lags, gamma_h, lag_counts

# Register the widget observers
register_observers([n_lags_entry, lag_resolution_entry, azimuth_range_entry], 
                    update_general_settings)

# General settings widget
def empirical_variogram_widget(settings_dict, plot_output_widget, style, layout):

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
    
    # Define the loading label
    loading_label = widgets.Label(value="")
    
    # Initialize the value dictionary
    general_settings = {
        "n_lags": n_lags_entry.value,
        "lag_resolution": lag_resolution_entry.value,
        "azimuth_range": azimuth_range_entry.value,
        "variable": empirical_variogram_dropdown.value,
    }

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
            model_name = VARIOGRAM_MODEL_WIDGET_MAP[model_selection]
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
                    default_value = get_variogram_defaults(survey, parameter)
                    # ---- Get the function argument name
                    arg_name = VARIOGRAM_MODEL_PARAMETER_MAP[parameter]["name"]
                    # ---- Get widget-type
                    widget_type = VARIOGRAM_MODEL_PARAMETER_MAP[parameter]["widget"]
                    # ---- Get step size
                    step_size = VARIOGRAM_MODEL_PARAMETER_MAP[parameter]["step"]
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
    copmute_theoretical_variogram_button.on_click(on_compute_button_clicked)

    # Vertically concatenate/bundle widgets
    widget_box = widgets.VBox([
        dropdown_variogram_model,
        variogram_model_widgets,
        copmute_theoretical_variogram_button
    ])

    # Return the widget and settings dictionary
    return widget_box, variogram_settings

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

def stitch_variogram_accordion(settings_widget, variogram_widget):

    # Initialize the accordion tabs
    variogram_accordion = widgets.Accordion(layout=dict(width="450px"))

    # Populate the accordion children
    variogram_accordion.children = [settings_widget, variogram_widget]#, optimization_widget]

    # Set the accordion tab titles
    variogram_accordion.set_title(0, "Empirical variogram", )
    variogram_accordion.set_title(1, "Theoretical variogram", )
    # variogram_accordion.set_title(2, "Optimize variogram parameters", )

    # Return the accordion
    return variogram_accordion

def on_compute_button_clicked(button):
    # parameters = variogram_widgets
    parameters = variogram_settings
    with plot_output_widget: 
        # out.clear_output(wait=True)  # Clear previous plot
        clear_output(wait=True)
        fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
        plot_theoretical_variogram(fig, ax, parameters)
        fig.canvas.toolbar_visible = False
        fig.canvas.header_visible = False
        fig.canvas.footer_visible = False
        fig.canvas.resizable = False
        plt.show(fig)
        # plt.close(fig)
        
def compute_theoretical_variogram(parameters):
    # Get the current input values
    arg_dict = {
        key: value.value if isinstance(value, (widgets.FloatText, widgets.Checkbox)) else value
        for key, value in parameters.items()
    }

    # pass into variogram
    
    distance_lags = lags
    return variogram(distance_lags, arg_dict)

def plot_theoretical_variogram(fig, ax, parameters):
    theoretical_variogram = compute_theoretical_variogram(parameters)
    
    ax.plot(lags, theoretical_variogram, linewidth = 3.0, c = "black")
    # plt.plot(lags, theoretical_variogram)
    fig.canvas.draw_idle()

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
# # ---- Theoretical variogram
# tab_contents_variogram, variogram_settings = theoretical_variogram_widget()
# # display(compute_empirical_variogram_button, plot_output_widget)

# # Initialize the accordion tabs
# # variogram_accordion = widgets.Accordion(layout=dict(width="550px"))
# # Create accordion
# variogram_accordion = stitch_variogram_accordion(tab_contents_settings, 
#                                                  tab_contents_variogram, )
# # Populate the accordion children
# # variogram_accordion.children = [tab_contents_settings]
# # Display the tabs
# # tab_contents_settings
# display(widgets.HBox([variogram_accordion, widgets.VBox([loading_label, plot_output_widget])]))


# def empirical_variogram_widget(settings_dict, plot_output_widget, style, layout):

#     # Define `matplotlib` backend
#     ipython_backend = get_ipython()
#     # ---- Set the backend (`%matplotlib widget`)
#     ipython_backend.run_line_magic("matplotlib", "widget")

#     # Create general settings value entry widgets
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
#                           description='Azimuth range', 
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
#     compute_empirical_variogram_button = widgets.Button(description="Compute empirical variogram",
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
#     # compute_empirical_variogram_button.on_click(click_empirical_variogram_button)

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
#                         update_general_settings)
    
#     # Vertically concatenate/bundle widgets
#     widget_box = widgets.VBox(
#         [n_lags_entry, lag_resolution_entry, azimuth_range_entry, empirical_variogram_dropdown,
#         compute_empirical_variogram_button]
#     )

#     # Return the widget and settings dictionary
#     return widget_box, general_settings

# widget_box