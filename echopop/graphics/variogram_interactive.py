import ipywidgets as widgets
import matplotlib.pyplot as plt
from typing import Dict
from IPython.display import display, clear_output

import pandas as pd
import numpy as np

from ..spatial.variogram import (
    get_variogram_arguments, empirical_variogram
)

from ..spatial.projection import transform_geometry

from ..spatial.transect import edit_transect_columns

########################
# Variogram model dropdown widget map
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

########################
# Empirical variogram computation function
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
    # Compute the empirical variogram
    lags, gamma_h, lag_counts, _ = compute_empirical_variogram()

    # Clear the previous plot
    variogram_plot_widget.clear_output()

    # Plot the new empirical variogram
    with variogram_plot_widget:
        fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
        plt.show(fig)

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

def register_observers(widgets, handler, names="value"):
    for widget in widgets:
        widget.observe(handler, names=names)

def empirical_variogram_widget(style, layout):

    # Initialize variogram settings for text widgets
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

def register_observers(widgets, handler, names='value'):
    for widget in widgets:
        widget.observe(handler, names=names)

# General settings widget
def general_settings_widget():

    # Create general settings value entry widgets
    # ---- Number of lags
    n_lags_entry = widgets.IntText(value=DEFAULT_SETTINGS_PARAMETERS["n_lags"],
                                   description="Number of lags",
                                   step=1, style=style, layout=layout)
    # ---- Lag (increment) resolution
    lag_resolution_entry =  widgets.FloatText(value=DEFAULT_SETTINGS_PARAMETERS["lag_resolution"], 
                                              description='Lag resolution', 
                                              step=DEFAULT_SETTINGS_PARAMETERS["lag_resolution"]/2, 
                                              style=style, layout=layout)
    # ---- Azimuth range 
    azimuth_range_entry = widgets.FloatText(value=DEFAULT_SETTINGS_PARAMETERS["azimuth_range"], 
                                            description='Azimuth Range', 
                                            step=DEFAULT_SETTINGS_PARAMETERS["azimuth_range"]/360, 
                                            style=style, layout=layout)
    
    # Add variable dropdown widget
    dropdown_variable = widgets.Dropdown(
        options=list(["Biomass (biomass density)", "Abundance (number density)"]),
        value="Biomass (biomass density)",
        description="Transect variable",
        disabled=False,
        style = {"description_width": "150px"}, 
        layout = {"width": "400px"}
    )    
    
    # Add the empirical variogram button
    compute_empirical_variogram_button = widgets.Button(description="Compute Empirical Variogram",
                                                        layout={"width": "405px"},
                                                        style={"font_weight": "bold"})
    
    # Initialize the value dictionary
    general_settings = {
        "n_lags": n_lags_entry.value,
        "lag_resolution": lag_resolution_entry.value,
        "azimuth_range": azimuth_range_entry.value,
        "variable": dropdown_variable.value,
    }

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
        [n_lags_entry, lag_resolution_entry, azimuth_range_entry, dropdown_variable,
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

    


def variogram_analysis_gui(transect_dict: Dict,
                           isobath_df: pd.DataFrame,
                           variogram_parameters: Dict,
                           settings_dict: Dict):

    # Define GUI widgets style and layout
    # ---- Style
    style = {"description_width": "150px"}
    # ---- Layout
    layout = {"width": "400px"}

    # Parameterize widgets
    # ---- Empirical variogram
    tab_contents_settings, general_settings = general_settings_widget()
    # ---- Theoretical variogram
    tab_contents_variogram, variogram_settings = theoretical_variogram_widget()
    # ---- Optimize variogram parameters
    placeholder_optimization = widgets.Label(value="Optimization Parameters tab will go here.")
    tab_contents_optimization = widgets.VBox([placeholder_optimization])

    # Update the settings dictionary
    settings_dict.update({
        "variable": (
            "biomass" if general_settings["variable"] == "Biomass (biomass density)" 
            else "abundance"
        ),
    })

    # Extract the transect data with the correct columns
    transect_data = edit_transect_columns(transect_dict, settings_dict)

    # Transform transect data geometry (standardized coordinates)
    transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

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

    # Initialize the Output plotting widget
    plot_output_widget = widgets.Output(layout=widgets.Layout(height="550px"))

    # Attach the plot to the Output widget
    with plot_output_widget:
        fig, ax = plot_empirical_variogram(lags, gamma_h, lag_counts)
        plt.show(fig)

    # Create accordion
    variogram_accordion = stitch_variogram_accordion(tab_contents_settings, 
                                                    tab_contents_variogram, 
                                                    tab_contents_optimization)

    # Display the tabs
    display(widgets.HBox([variogram_accordion, plot_output_widget]))






# import tkinter as tk
# from tkinter import ttk
# from tkinter import messagebox
# from tkinter import simpledialog
# from tkinter import filedialog
# from typing import Type  # you have to import Type
# import ipywidgets as widgets
# from IPython.display import display, clear_output
# import matplotlib.pyplot as plt
# from echopop.spatial.variogram import get_variogram_arguments
# from echopop.analysis import (
#     variogram_analysis,
# )


# from matplotlib.figure import Figure
# from matplotlib.backends.backend_agg import FigureCanvasAgg

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
#     # plt.close(fig)

#     return fig
    
#     # # Create a list of annotated points
#     # annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
#     #                     bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
#     #                     arrowprops=dict(arrowstyle="->", color="black", lw=1, linestyle="-"))
#     # annot.set_visible(False)

#     # def update_annot(ind):
#     #     pos = scatter.get_offsets()[ind["ind"][0]]
#     #     annot.xy = pos
#     #     text = f"Lag={pos[0]:.2f}\nγ={pos[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}"
#     #     annot.set_text(text)
#     #     annot.get_bbox_patch().set_alpha(1.0)

#     # def on_hover(event):
#     #     vis = annot.get_visible()
#     #     if event.inaxes == ax:
#     #         cont, ind = scatter.contains(event)
#     #         if cont:
#     #             update_annot(ind)
#     #             annot.set_visible(True)
#     #             fig.canvas.draw_idle()
#     #         else:
#     #             if vis:
#     #                 annot.set_visible(False)
#     #                 fig.canvas.draw_idle()

#     # fig.canvas.mpl_connect("motion_notify_event", on_hover)

#     # return fig


# style = {'description_width': '150px'}
# layout = {'width': '400px'}
# # Dropdown widget for variogram model selection
# dropdown_variogram_model = widgets.Dropdown(
#     options=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()),
#     value="Bessel-exponential",
#     description='Variogram Model',
#     disabled=False,
#     style=style, layout=layout
# )
# # Output widget to display parameter fields dynamically
# output_variogram_params = widgets.Output(layout={'height': 'auto', "width": "auto"})
# output_optimization_params = widgets.Output()

# # Function to update variogram parameters based on dropdown selection
# def update_variogram_params(change):
#     with output_variogram_params:
#         clear_output(wait=True)
#         selected_model = change['new']
#         print(f"Selected Model Type: {type(selected_model)}")
#         print(f"Selected Model Value: {selected_model}")

#         model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
#         print(f"{model_names}")

#         function_arguments, _ = get_variogram_arguments(model_names)
#         # ---- convert to a dictionary
#         args = dict(function_arguments)

#         print(f"{args.keys()}")

#         # Clear existing widgets (if any)
#         output_variogram_params.clear_output()

#         # Display widgets for each parameter
#         for param in args.values():
#             arg = param.name

#             if arg != "distance_lags":  # Adjust as needed
#                 default_value = get_variogram_defaults(survey, arg)
#                 # ---- get the name from the key
#                 arg_name = param_key[arg]["name"]
#                 widget_type = param_key[arg]["widget"]

#                 if widget_type == "checkbox":
#                     checkbox = widgets.Checkbox(value=default_value, description=arg_name, style={'description_width': 'initial'})
#                     display(checkbox)
#                 # Example: Create widgets (Checkbox and FloatText) for each parameter
#                 elif widget_type == "entry":
#                     entry = widgets.FloatText(value=default_value, description=arg_name, style=style, layout=layout)
#                     display(entry)

# # Attach the observer function to the dropdown
# dropdown_variogram_model.observe(update_variogram_params, names='value')
# update_variogram_params({'new': dropdown_variogram_model.value})
# placeholder_optimization = widgets.Label(value='Optimization Parameters tab will go here.')

# tab_contents_variogram = widgets.VBox([dropdown_variogram_model, output_variogram_params])
# tab_contents_optimization = widgets.VBox([placeholder_optimization])

# accordion_layout = {'width': '450px'}
# tab = widgets.Accordion(layout=accordion_layout)
# tab.children = [tab_contents_variogram, tab_contents_optimization]

# # Set titles for tabs using Markdown for better formatting
# tab.set_title(0, 'Variogram Parameters', )
# tab.set_title(1, 'Optimization Parameters')
# # Output widget for the plot
# plot_output = widgets.Output()

# layout = widgets.HBox([tab, plot_output])

# # Display the layout
# display(layout)

# # Function to display the plot in plot_output
# def display_plot():
#     plot_output.clear_output(wait=True)
#     with plot_output:
#         fig = plot_variogram()
#         display(fig.canvas)

# # Initial plot rendering
# display_plot()
# # fig, ax = plot_variogram()
# # Display the tab widget
# # plot_output = widgets.Output()
# # with plot_output:
# #     display(fig.canvas)

# # w = widgets.interactive(plot_variogram, tab)

# # left_panel = widgets.VBox([tab])
# # plot_output
# # widgets.HBox([right_panel])

# # widgets.HBox([tab, plot_output])


# # ---- Get `transect_data`
# lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
#     transect_data, variogram_parameters, settings_dict
# )

# survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
#                 survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")

# def get_variogram_defaults(survey, argument):
#     DEFAULT_PARAMETERS = {
#         "nugget": 0.0,
#         "sill": 1.0,
#         "hole_effect_range": 0.0,
#         "decay_power": 1.5,
#         "enhance_semivariance": False,
#     }
#     # ---- get variogram config
#     if "model_config" in survey.input["statistics"]["variogram"]:
#         if argument in survey.input["statistics"]["variogram"]["model_config"]:
#             return survey.input["statistics"]["variogram"]["model_config"][argument]
#         else:
#             return DEFAULT_PARAMETERS[argument]

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
# selected_model = "Bessel-exponential"
# param_key = {
#     "nugget": dict(name="Nugget", widget="entry"),
#     "sill": dict(name="Sill", widget="entry"),
#     "correlation_range": dict(name="Correlation range", widget="entry"),
#     "decay_power": dict(name="Decay power", widget="entry"),
#     "hole_effect_range": dict(name="Hole effect range", widget="entry"),
#     "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox")
# }

# dropdown_variogram_model = widgets.Dropdown(
#     options=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()),
#     value="Bessel-exponential",
#     description='Variogram Model:',
#     disabled=False,
# )
# dropdown_variogram_model.value

# def update_variogram_params(change):
#     with output_variogram_params:
#         clear_output(wait=True)
#         selected_model = change.new
#         print(f"Selected Model Type: {type(selected_model)}")
#         print(f"Selected Model Value: {selected_model}")

#         model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
#         print(f"{model_names}")

#         function_arguments, _ = get_variogram_arguments(model_names)
#         # ---- convert to a dictionary
#         args = dict(function_arguments)

#         print(f"{args.keys()}")

#         # Clear existing widgets (if any)
#         output_variogram_params.clear_output()

#         # Display widgets for each parameter
#         for param in args.values():
#             arg = param.name

#             if arg != "distance_lags":  # Adjust as needed
#                 default_value = get_variogram_defaults(survey, arg)
#                 # ---- get the name from the key
#                 arg_name = param_key[arg]["name"]
#                 widget_type = param_key[arg]["widget"]

#                 if widget_type == "checkbox":
#                     checkbox = widgets.Checkbox(value=default_value, description=arg_name)
#                     display(checkbox)
#                 # Example: Create widgets (Checkbox and FloatText) for each parameter
#                 elif widget_type == "entry":
#                     entry = widgets.FloatText(value=default_value, description=arg_name)
#                     display(entry)
# model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
# function_arguments, _ = get_variogram_arguments(model_names)
# args = dict(function_arguments)
# param = args["nugget"]

# for param in args.values():
#     # Assuming param_name is the parameter name and param_description is its description
#     # Example: Create widgets (Checkbox and Entry) for each parameter
#     arg = param.name
#     argument = arg
#     if arg != "distance_lags":
#         default_value = get_variogram_defaults(survey, argument)
#         # ---- get the name from the key
#         arg_name = param_key[arg]
#         print(arg_name)

#     # ---- get the name from the key
#     arg_name = param_key[arg]

#     if arg != "distance_lags":
#         print(arg_name)

# # selected_model = change.new
# with output_variogram_params:  
#     model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
#     function_arguments, _ = get_variogram_arguments(model_names)
#     # ---- convert to a dictionary
#     args = dict(function_arguments)
    
#     for param in args.values():
#         # Assuming param_name is the parameter name and param_description is its description
#         # Example: Create widgets (Checkbox and Entry) for each parameter
#         arg = param.name
#         # ---- get the name from the key
#         arg_name = param_key[arg]
#         checkbox = widgets.Checkbox(value=True, description=arg_name)
#         display(checkbox)
#         if param_name == "decay_power":  # Example condition for specific handling
#             entry = widgets.FloatText(description="Value:")
#             display(entry)
# VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
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

# # Create the main application window
# root = tk.Tk()
# root.title("Fit Variogram GUI")

# # Set minimum size of the window
# root.minsize(300, 100)  # Set minimum width and height

# # Create tabs using ttk.Notebook
# notebook = ttk.Notebook(root)
# notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)


# # Add input fields for variogram parameters
# variogram_fields = [
#     ("Number of Lags", "n_lags"),
#     ("Lag Resolution", "lag_resolution"),
#     ("Max Range", "max_range"),
#     ("Sill", "sill"),
#     ("Nugget", "nugget"),
#     ("Hole Effect Range", "hole_effect_range"),
#     ("Correlation Range", "correlation_range"),
#     ("Enhance Semivariance", "enhance_semivariance"),
#     ("Decay Power", "decay_power")
# ]

# # Add input fields for optimization parameters
# optimization_fields = [
#     ("Max Function Evaluations", "max_fun_evaluations"),
#     ("Cost Function Tolerance", "cost_fun_tolerance"),
#     ("Solution Tolerance", "solution_tolerance"),
#     ("Gradient Tolerance", "gradient_tolerance"),
#     ("Finite Step Size", "finite_step_size"),
#     ("Trust Region Solver", "trust_region_solver"),
#     ("X Scale", "x_scale"),
#     ("Jacobian Approx", "jacobian_approx")
# ]

# def on_selection(event):
#     """
#     """
#     selected_model = combobox.get()
#     get_variogram_model(selected_model)

# def get_variogram_model(selected_model):
#     """
#     """

#     # Get the corresponding value from the `WIDGET_MAP`
#     model = WIDGET_MAP.get(selected_model)

#     # Begin processing
#     if isinstance(model, list):
#         print(f"Processing composite model: {selected_model}")
#     else:
#         print(f"Processing model: {selected_model}")

# def generate_checkboxes(parent, fit_params):
#     for i, param in enumerate(fit_params):
#         var = tk.IntVar()
#         checkbox = ttk.Checkbutton(parent, text=param, variable=var)
#         checkbox.grid(row=i, column=0, padx=10, pady=5, sticky=tk.W)

# # Function to update fit parameters based on selected model
# def update_fit_parameters(event):
#     selected_model = combobox.get()
#     print(f"Selected Model: {selected_model}")

#     # Clear existing widgets in the fit parameters tab
#     for widget in variogram_frame.winfo_children():
#         widget.destroy()

#     # Generate fit parameters based on selected model
#     if selected_model in WIDGET_MAP:
#         fit_params = WIDGET_MAP[selected_model]
#         for field, var_name in variogram_fields:
#             label = ttk.Label(variogram_frame, text=field)
#             label.pack(padx=10, pady=2)
#             if var_name == "enhance_semivariance":
#                 var = tk.BooleanVar(value=True)
#                 checkbutton = ttk.Checkbutton(variogram_frame, variable=var)
#                 checkbutton.pack(padx=10, pady=2)
#             else:
#                 entry = ttk.Entry(variogram_frame)
#                 entry.pack(padx=10, pady=2)

#     else:
#         print(f"No fit parameters defined for model: {selected_model}")


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


# # Create the main application window
# root = tk.Tk()
# root.title("Fit Variogram GUI")

# # Set minimum size of the window
# root.minsize(300, 100)  # Set minimum width and height

# # Create tabs using ttk.Notebook
# notebook = ttk.Notebook(root)
# notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# # Create frames for each tab
# variogram_frame = ttk.Frame(notebook)
# optimization_frame = ttk.Frame(notebook)

# # Add tabs to the notebook
# notebook.add(variogram_frame, text='Variogram Parameters')
# notebook.add(optimization_frame, text='Optimization Parameters')

# # Dropdown widget
# combobox = ttk.Combobox(root, values=list(WIDGET_MAP.keys()), state = "readonly", width=20)
# combobox.bind("<<ComboboxSelected>>", on_selection)
# # Bind selection event to update fit parameters
# combobox.bind("<<ComboboxSelected>>", update_fit_parameters)
# combobox.pack(padx=10, pady=5)

# for field, var_name in variogram_fields:
#     label = ttk.Label(variogram_frame, text=field)
#     label.pack(padx=10, pady=2)
#     if var_name in ["enhance_semivariance"]:
#         var = tk.BooleanVar(value=True)
#         checkbutton = ttk.Checkbutton(variogram_frame, variable=var)
#         checkbutton.pack(padx=10, pady=2)
#     else:
#         entry = ttk.Entry(variogram_frame)
#         entry.pack(padx=10, pady=2)

# # Initialize optimization frame
# create_optimization_frame()

# # Add a button to call the fit_variogram function
# button = ttk.Button(root, text="Fit Variogram")
# button.pack(padx=10, pady=10)

# # Start the main event loop
# root.mainloop()


# # Example WIDGET_MAP and fit_parameters setup
# WIDGET_MAP = {
#     "Bessel-exponential": ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"],
#     "Bessel": ["nugget", "sill", "correlation_range"],
#     "Exponential": ["nugget", "sill"],
#     "Gaussian": ["nugget", "sill", "correlation_range"],
#     "Linear": ["sill"],
#     "Sinc": ["sill", "correlation_range"],
#     "Spherical": ["nugget", "sill", "correlation_range", "hole_effect_range"],
#     "Bessel-Gaussian": ["nugget", "sill", "correlation_range"],
#     "Cosine-exponential": ["nugget", "sill", "decay_power"],
#     "Cosine-Gaussian": ["nugget", "sill", "correlation_range"],
#     "Exponential-linear": ["sill", "decay_power"],
#     "Gaussian-linear": ["sill", "correlation_range"]
# }

# import tkinter as tk
# from tkinter import ttk

# # Example WIDGET_MAP and fit_parameters setup
# WIDGET_MAP = {
#     "Bessel-exponential": ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"],
#     "Bessel": ["nugget", "sill", "correlation_range"],
#     "Exponential": ["nugget", "sill"],
#     "Gaussian": ["nugget", "sill", "correlation_range"],
#     "Linear": ["sill"],
#     "Sinc": ["sill", "correlation_range"],
#     "Spherical": ["nugget", "sill", "correlation_range", "hole_effect_range"],
#     "Bessel-Gaussian": ["nugget", "sill", "correlation_range"],
#     "Cosine-exponential": ["nugget", "sill", "decay_power"],
#     "Cosine-Gaussian": ["nugget", "sill", "correlation_range"],
#     "Exponential-linear": ["sill", "decay_power"],
#     "Gaussian-linear": ["sill", "correlation_range"]
# }

# # Create the main application window
# root = tk.Tk()
# root.title("Fit Variogram GUI")

# # Set minimum size of the window
# root.minsize(400, 300)  # Adjust as per your minimum size requirement

# # Create tabs using ttk.Notebook
# notebook = ttk.Notebook(root)
# notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# # Create frames for each tab
# variogram_frame = ttk.Frame(notebook)
# optimization_frame = ttk.Frame(notebook)

# # Add tabs to the notebook
# notebook.add(variogram_frame, text='Variogram Parameters')
# notebook.add(optimization_frame, text='Optimization Parameters')

# # Dropdown widget
# combobox = ttk.Combobox(root, values=list(WIDGET_MAP.keys()), state="readonly", width=20)
# combobox.pack(padx=10, pady=5)

# # Function to update fit parameters based on selected model
# def update_fit_parameters(event):
#     selected_model = combobox.get()
#     print(f"Selected Model: {selected_model}")

#     # Clear existing widgets in the variogram frame
#     for widget in variogram_frame.winfo_children():
#         widget.destroy()

#     # Generate checkboxes and entry fields based on selected model
#     if selected_model in WIDGET_MAP:
#         fit_params = WIDGET_MAP[selected_model]
#         for param in fit_params:
#             var = tk.BooleanVar(value=True)
#             checkbox = ttk.Checkbutton(variogram_frame, text=param, variable=var)
#             checkbox.pack(side=tk.LEFT, padx=5, pady=5)

#             # Add an entry field next to each checkbox
#             entry = ttk.Entry(variogram_frame, width=10)
#             entry.pack(side=tk.LEFT, padx=5, pady=5)
#     else:
#         print(f"No fit parameters defined for model: {selected_model}")

# # Bind selection event to update fit parameters
# combobox.bind("<<ComboboxSelected>>", update_fit_parameters)

# # Function to create the optimization parameters frame
# def create_optimization_frame():
#     optimization_fields = [
#         ("Max Fun Evaluations", "max_fun_evaluations"),
#         ("Cost Fun Tolerance", "cost_fun_tolerance"),
#         ("Solution Tolerance", "solution_tolerance"),
#         ("Gradient Tolerance", "gradient_tolerance"),
#         ("Finite Step Size", "finite_step_size"),
#         ("Trust Region Solver", "trust_region_solver"),
#         ("X Scale", "x_scale"),
#         ("Jacobian Approx", "jacobian_approx"),
#         ("Force Lag Zero", "force_lag_zero"),
#         ("Verbose", "verbose")
#     ]
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

# # Initialize optimization frame
# create_optimization_frame()

# # Add a button to call the fit_variogram function
# button = ttk.Button(root, text="Fit Variogram")
# button.pack(padx=10, pady=10)

# # Start the main event loop
# root.mainloop()
