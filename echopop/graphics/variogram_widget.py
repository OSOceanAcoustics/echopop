import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import simpledialog
from tkinter import filedialog
from typing import Type  # you have to import Type
import ipywidgets as widgets
from IPython.display import display, clear_output
import matplotlib.pyplot as plt
from echopop.spatial.variogram import get_variogram_arguments
from echopop.analysis import (
    variogram_analysis,
)


from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

def plot_variogram():
    s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
    x = lags; y = gamma_h

    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(x, y, s=s(lag_counts))
    ax.set_xlabel('Lag Distance')
    ax.set_ylabel('Semivariance')
    ax.set_title('Empirical Variogram')
    ax.grid(True)
    plt.tight_layout()
    # plt.close(fig)

    return fig
    
    # # Create a list of annotated points
    # annot = ax.annotate("", xy=(0, 0), xytext=(10, 10), textcoords="offset points",
    #                     bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
    #                     arrowprops=dict(arrowstyle="->", color="black", lw=1, linestyle="-"))
    # annot.set_visible(False)

    # def update_annot(ind):
    #     pos = scatter.get_offsets()[ind["ind"][0]]
    #     annot.xy = pos
    #     text = f"Lag={pos[0]:.2f}\nγ={pos[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}"
    #     annot.set_text(text)
    #     annot.get_bbox_patch().set_alpha(1.0)

    # def on_hover(event):
    #     vis = annot.get_visible()
    #     if event.inaxes == ax:
    #         cont, ind = scatter.contains(event)
    #         if cont:
    #             update_annot(ind)
    #             annot.set_visible(True)
    #             fig.canvas.draw_idle()
    #         else:
    #             if vis:
    #                 annot.set_visible(False)
    #                 fig.canvas.draw_idle()

    # fig.canvas.mpl_connect("motion_notify_event", on_hover)

    # return fig


style = {'description_width': '150px'}
layout = {'width': '400px'}
# Dropdown widget for variogram model selection
dropdown_variogram_model = widgets.Dropdown(
    options=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()),
    value="Bessel-exponential",
    description='Variogram Model',
    disabled=False,
    style=style, layout=layout
)
# Output widget to display parameter fields dynamically
output_variogram_params = widgets.Output(layout={'height': 'auto', "width": "auto"})
output_optimization_params = widgets.Output()

# Function to update variogram parameters based on dropdown selection
def update_variogram_params(change):
    with output_variogram_params:
        clear_output(wait=True)
        selected_model = change['new']
        print(f"Selected Model Type: {type(selected_model)}")
        print(f"Selected Model Value: {selected_model}")

        model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
        print(f"{model_names}")

        function_arguments, _ = get_variogram_arguments(model_names)
        # ---- convert to a dictionary
        args = dict(function_arguments)

        print(f"{args.keys()}")

        # Clear existing widgets (if any)
        output_variogram_params.clear_output()

        # Display widgets for each parameter
        for param in args.values():
            arg = param.name

            if arg != "distance_lags":  # Adjust as needed
                default_value = get_variogram_defaults(survey, arg)
                # ---- get the name from the key
                arg_name = param_key[arg]["name"]
                widget_type = param_key[arg]["widget"]

                if widget_type == "checkbox":
                    checkbox = widgets.Checkbox(value=default_value, description=arg_name, style={'description_width': 'initial'})
                    display(checkbox)
                # Example: Create widgets (Checkbox and FloatText) for each parameter
                elif widget_type == "entry":
                    entry = widgets.FloatText(value=default_value, description=arg_name, style=style, layout=layout)
                    display(entry)

# Attach the observer function to the dropdown
dropdown_variogram_model.observe(update_variogram_params, names='value')
update_variogram_params({'new': dropdown_variogram_model.value})
placeholder_optimization = widgets.Label(value='Optimization Parameters tab will go here.')

tab_contents_variogram = widgets.VBox([dropdown_variogram_model, output_variogram_params])
tab_contents_optimization = widgets.VBox([placeholder_optimization])

accordion_layout = {'width': '450px'}
tab = widgets.Accordion(layout=accordion_layout)
tab.children = [tab_contents_variogram, tab_contents_optimization]

# Set titles for tabs using Markdown for better formatting
tab.set_title(0, 'Variogram Parameters', )
tab.set_title(1, 'Optimization Parameters')
# Output widget for the plot
plot_output = widgets.Output()

layout = widgets.HBox([tab, plot_output])

# Display the layout
display(layout)

# Function to display the plot in plot_output
def display_plot():
    plot_output.clear_output(wait=True)
    with plot_output:
        fig = plot_variogram()
        display(fig.canvas)

# Initial plot rendering
display_plot()
# fig, ax = plot_variogram()
# Display the tab widget
# plot_output = widgets.Output()
# with plot_output:
#     display(fig.canvas)

# w = widgets.interactive(plot_variogram, tab)

# left_panel = widgets.VBox([tab])
# plot_output
# widgets.HBox([right_panel])

# widgets.HBox([tab, plot_output])


# ---- Get `transect_data`
lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
    transect_data, variogram_parameters, settings_dict
)

survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")

def get_variogram_defaults(survey, argument):
    DEFAULT_PARAMETERS = {
        "nugget": 0.0,
        "sill": 1.0,
        "hole_effect_range": 0.0,
        "decay_power": 1.5,
        "enhance_semivariance": False,
    }
    # ---- get variogram config
    if "model_config" in survey.input["statistics"]["variogram"]:
        if argument in survey.input["statistics"]["variogram"]["model_config"]:
            return survey.input["statistics"]["variogram"]["model_config"][argument]
        else:
            return DEFAULT_PARAMETERS[argument]

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
selected_model = "Bessel-exponential"
param_key = {
    "nugget": dict(name="Nugget", widget="entry"),
    "sill": dict(name="Sill", widget="entry"),
    "correlation_range": dict(name="Correlation range", widget="entry"),
    "decay_power": dict(name="Decay power", widget="entry"),
    "hole_effect_range": dict(name="Hole effect range", widget="entry"),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox")
}

dropdown_variogram_model = widgets.Dropdown(
    options=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()),
    value="Bessel-exponential",
    description='Variogram Model:',
    disabled=False,
)
dropdown_variogram_model.value

def update_variogram_params(change):
    with output_variogram_params:
        clear_output(wait=True)
        selected_model = change.new
        print(f"Selected Model Type: {type(selected_model)}")
        print(f"Selected Model Value: {selected_model}")

        model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
        print(f"{model_names}")

        function_arguments, _ = get_variogram_arguments(model_names)
        # ---- convert to a dictionary
        args = dict(function_arguments)

        print(f"{args.keys()}")

        # Clear existing widgets (if any)
        output_variogram_params.clear_output()

        # Display widgets for each parameter
        for param in args.values():
            arg = param.name

            if arg != "distance_lags":  # Adjust as needed
                default_value = get_variogram_defaults(survey, arg)
                # ---- get the name from the key
                arg_name = param_key[arg]["name"]
                widget_type = param_key[arg]["widget"]

                if widget_type == "checkbox":
                    checkbox = widgets.Checkbox(value=default_value, description=arg_name)
                    display(checkbox)
                # Example: Create widgets (Checkbox and FloatText) for each parameter
                elif widget_type == "entry":
                    entry = widgets.FloatText(value=default_value, description=arg_name)
                    display(entry)
model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
function_arguments, _ = get_variogram_arguments(model_names)
args = dict(function_arguments)
param = args["nugget"]

for param in args.values():
    # Assuming param_name is the parameter name and param_description is its description
    # Example: Create widgets (Checkbox and Entry) for each parameter
    arg = param.name
    argument = arg
    if arg != "distance_lags":
        default_value = get_variogram_defaults(survey, argument)
        # ---- get the name from the key
        arg_name = param_key[arg]
        print(arg_name)

    # ---- get the name from the key
    arg_name = param_key[arg]

    if arg != "distance_lags":
        print(arg_name)

# selected_model = change.new
with output_variogram_params:  
    model_names = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
    function_arguments, _ = get_variogram_arguments(model_names)
    # ---- convert to a dictionary
    args = dict(function_arguments)
    
    for param in args.values():
        # Assuming param_name is the parameter name and param_description is its description
        # Example: Create widgets (Checkbox and Entry) for each parameter
        arg = param.name
        # ---- get the name from the key
        arg_name = param_key[arg]
        checkbox = widgets.Checkbox(value=True, description=arg_name)
        display(checkbox)
        if param_name == "decay_power":  # Example condition for specific handling
            entry = widgets.FloatText(description="Value:")
            display(entry)
VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
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

optimization_fields = [
    ("Maximum function evaluations", "max_fun_evaluations"),
    ("Cost function tolerance", "cost_fun_tolerance"),
    ("Solution tolerance", "solution_tolerance"),
    ("Gradient tolerance", "gradient_tolerance"),
    ("Finite differences step size", "finite_step_size"),
    ("Trust Region solver method", "trust_region_solver"),
    ("Characteristic x-scale", "x_scale"),
    ("Jacobian approximation method", "jacobian_approx")
]

# Create the main application window
root = tk.Tk()
root.title("Fit Variogram GUI")

# Set minimum size of the window
root.minsize(300, 100)  # Set minimum width and height

# Create tabs using ttk.Notebook
notebook = ttk.Notebook(root)
notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)


# Add input fields for variogram parameters
variogram_fields = [
    ("Number of Lags", "n_lags"),
    ("Lag Resolution", "lag_resolution"),
    ("Max Range", "max_range"),
    ("Sill", "sill"),
    ("Nugget", "nugget"),
    ("Hole Effect Range", "hole_effect_range"),
    ("Correlation Range", "correlation_range"),
    ("Enhance Semivariance", "enhance_semivariance"),
    ("Decay Power", "decay_power")
]

# Add input fields for optimization parameters
optimization_fields = [
    ("Max Function Evaluations", "max_fun_evaluations"),
    ("Cost Function Tolerance", "cost_fun_tolerance"),
    ("Solution Tolerance", "solution_tolerance"),
    ("Gradient Tolerance", "gradient_tolerance"),
    ("Finite Step Size", "finite_step_size"),
    ("Trust Region Solver", "trust_region_solver"),
    ("X Scale", "x_scale"),
    ("Jacobian Approx", "jacobian_approx")
]

def on_selection(event):
    """
    """
    selected_model = combobox.get()
    get_variogram_model(selected_model)

def get_variogram_model(selected_model):
    """
    """

    # Get the corresponding value from the `WIDGET_MAP`
    model = WIDGET_MAP.get(selected_model)

    # Begin processing
    if isinstance(model, list):
        print(f"Processing composite model: {selected_model}")
    else:
        print(f"Processing model: {selected_model}")

def generate_checkboxes(parent, fit_params):
    for i, param in enumerate(fit_params):
        var = tk.IntVar()
        checkbox = ttk.Checkbutton(parent, text=param, variable=var)
        checkbox.grid(row=i, column=0, padx=10, pady=5, sticky=tk.W)

# Function to update fit parameters based on selected model
def update_fit_parameters(event):
    selected_model = combobox.get()
    print(f"Selected Model: {selected_model}")

    # Clear existing widgets in the fit parameters tab
    for widget in variogram_frame.winfo_children():
        widget.destroy()

    # Generate fit parameters based on selected model
    if selected_model in WIDGET_MAP:
        fit_params = WIDGET_MAP[selected_model]
        for field, var_name in variogram_fields:
            label = ttk.Label(variogram_frame, text=field)
            label.pack(padx=10, pady=2)
            if var_name == "enhance_semivariance":
                var = tk.BooleanVar(value=True)
                checkbutton = ttk.Checkbutton(variogram_frame, variable=var)
                checkbutton.pack(padx=10, pady=2)
            else:
                entry = ttk.Entry(variogram_frame)
                entry.pack(padx=10, pady=2)

    else:
        print(f"No fit parameters defined for model: {selected_model}")


# Function to create the optimization parameters frame
def create_optimization_frame():
    for field, var_name in optimization_fields:
        label = ttk.Label(optimization_frame, text=field)
        label.pack(padx=10, pady=2)
        if var_name in ["force_lag_zero", "verbose"]:
            var = tk.BooleanVar(value=True)
            checkbutton = ttk.Checkbutton(optimization_frame, variable=var)
            checkbutton.pack(padx=10, pady=2)
        else:
            entry = ttk.Entry(optimization_frame)
            entry.pack(padx=10, pady=2)


# Create the main application window
root = tk.Tk()
root.title("Fit Variogram GUI")

# Set minimum size of the window
root.minsize(300, 100)  # Set minimum width and height

# Create tabs using ttk.Notebook
notebook = ttk.Notebook(root)
notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Create frames for each tab
variogram_frame = ttk.Frame(notebook)
optimization_frame = ttk.Frame(notebook)

# Add tabs to the notebook
notebook.add(variogram_frame, text='Variogram Parameters')
notebook.add(optimization_frame, text='Optimization Parameters')

# Dropdown widget
combobox = ttk.Combobox(root, values=list(WIDGET_MAP.keys()), state = "readonly", width=20)
combobox.bind("<<ComboboxSelected>>", on_selection)
# Bind selection event to update fit parameters
combobox.bind("<<ComboboxSelected>>", update_fit_parameters)
combobox.pack(padx=10, pady=5)

for field, var_name in variogram_fields:
    label = ttk.Label(variogram_frame, text=field)
    label.pack(padx=10, pady=2)
    if var_name in ["enhance_semivariance"]:
        var = tk.BooleanVar(value=True)
        checkbutton = ttk.Checkbutton(variogram_frame, variable=var)
        checkbutton.pack(padx=10, pady=2)
    else:
        entry = ttk.Entry(variogram_frame)
        entry.pack(padx=10, pady=2)

# Initialize optimization frame
create_optimization_frame()

# Add a button to call the fit_variogram function
button = ttk.Button(root, text="Fit Variogram")
button.pack(padx=10, pady=10)

# Start the main event loop
root.mainloop()


# Example WIDGET_MAP and fit_parameters setup
WIDGET_MAP = {
    "Bessel-exponential": ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"],
    "Bessel": ["nugget", "sill", "correlation_range"],
    "Exponential": ["nugget", "sill"],
    "Gaussian": ["nugget", "sill", "correlation_range"],
    "Linear": ["sill"],
    "Sinc": ["sill", "correlation_range"],
    "Spherical": ["nugget", "sill", "correlation_range", "hole_effect_range"],
    "Bessel-Gaussian": ["nugget", "sill", "correlation_range"],
    "Cosine-exponential": ["nugget", "sill", "decay_power"],
    "Cosine-Gaussian": ["nugget", "sill", "correlation_range"],
    "Exponential-linear": ["sill", "decay_power"],
    "Gaussian-linear": ["sill", "correlation_range"]
}

import tkinter as tk
from tkinter import ttk

# Example WIDGET_MAP and fit_parameters setup
WIDGET_MAP = {
    "Bessel-exponential": ["nugget", "sill", "correlation_range", "hole_effect_range", "decay_power"],
    "Bessel": ["nugget", "sill", "correlation_range"],
    "Exponential": ["nugget", "sill"],
    "Gaussian": ["nugget", "sill", "correlation_range"],
    "Linear": ["sill"],
    "Sinc": ["sill", "correlation_range"],
    "Spherical": ["nugget", "sill", "correlation_range", "hole_effect_range"],
    "Bessel-Gaussian": ["nugget", "sill", "correlation_range"],
    "Cosine-exponential": ["nugget", "sill", "decay_power"],
    "Cosine-Gaussian": ["nugget", "sill", "correlation_range"],
    "Exponential-linear": ["sill", "decay_power"],
    "Gaussian-linear": ["sill", "correlation_range"]
}

# Create the main application window
root = tk.Tk()
root.title("Fit Variogram GUI")

# Set minimum size of the window
root.minsize(400, 300)  # Adjust as per your minimum size requirement

# Create tabs using ttk.Notebook
notebook = ttk.Notebook(root)
notebook.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

# Create frames for each tab
variogram_frame = ttk.Frame(notebook)
optimization_frame = ttk.Frame(notebook)

# Add tabs to the notebook
notebook.add(variogram_frame, text='Variogram Parameters')
notebook.add(optimization_frame, text='Optimization Parameters')

# Dropdown widget
combobox = ttk.Combobox(root, values=list(WIDGET_MAP.keys()), state="readonly", width=20)
combobox.pack(padx=10, pady=5)

# Function to update fit parameters based on selected model
def update_fit_parameters(event):
    selected_model = combobox.get()
    print(f"Selected Model: {selected_model}")

    # Clear existing widgets in the variogram frame
    for widget in variogram_frame.winfo_children():
        widget.destroy()

    # Generate checkboxes and entry fields based on selected model
    if selected_model in WIDGET_MAP:
        fit_params = WIDGET_MAP[selected_model]
        for param in fit_params:
            var = tk.BooleanVar(value=True)
            checkbox = ttk.Checkbutton(variogram_frame, text=param, variable=var)
            checkbox.pack(side=tk.LEFT, padx=5, pady=5)

            # Add an entry field next to each checkbox
            entry = ttk.Entry(variogram_frame, width=10)
            entry.pack(side=tk.LEFT, padx=5, pady=5)
    else:
        print(f"No fit parameters defined for model: {selected_model}")

# Bind selection event to update fit parameters
combobox.bind("<<ComboboxSelected>>", update_fit_parameters)

# Function to create the optimization parameters frame
def create_optimization_frame():
    optimization_fields = [
        ("Max Fun Evaluations", "max_fun_evaluations"),
        ("Cost Fun Tolerance", "cost_fun_tolerance"),
        ("Solution Tolerance", "solution_tolerance"),
        ("Gradient Tolerance", "gradient_tolerance"),
        ("Finite Step Size", "finite_step_size"),
        ("Trust Region Solver", "trust_region_solver"),
        ("X Scale", "x_scale"),
        ("Jacobian Approx", "jacobian_approx"),
        ("Force Lag Zero", "force_lag_zero"),
        ("Verbose", "verbose")
    ]
    for field, var_name in optimization_fields:
        label = ttk.Label(optimization_frame, text=field)
        label.pack(padx=10, pady=2)
        if var_name in ["force_lag_zero", "verbose"]:
            var = tk.BooleanVar(value=True)
            checkbutton = ttk.Checkbutton(optimization_frame, variable=var)
            checkbutton.pack(padx=10, pady=2)
        else:
            entry = ttk.Entry(optimization_frame)
            entry.pack(padx=10, pady=2)

# Initialize optimization frame
create_optimization_frame()

# Add a button to call the fit_variogram function
button = ttk.Button(root, text="Fit Variogram")
button.pack(padx=10, pady=10)

# Start the main event loop
root.mainloop()
