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


transect_dict =  self.analysis["transect"]
settings_dict =  self.analysis["settings"]["variogram"]
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]

# Extract specific variogram parameters
# ---- Number of lags
n_lags = variogram_parameters["n_lags"]
# ---- Lag resolution
lag_resolution = variogram_parameters["lag_resolution"]

# Compute the lag distances
distance_lags = np.arange(1, n_lags) * lag_resolution
# ---- Add to the `variogram_parameters` dictionary
variogram_parameters["distance_lags"] = distance_lags
# ---- Update the max range parameter, if necessary
max_range = lag_resolution * n_lags


optimization_parameters["solution_tolerance"] = 1e-4
optimization_parameters["cost_fun_tolerance"] = 1e-6
optimization_parameters["gradient_tolerance"] = 1e-4
optimization_parameters["jacobian_approx"] = "central"
# Generate the optimization settings dictionary
optimization_settings = create_optimization_options(
    fit_parameters,
    variogram_parameters,
    initial_values=initial_values,
    lower_bounds=lower_bounds,
    upper_bounds=upper_bounds,
    model=variogram_parameters["model"],
    **optimization_parameters,
)

# Validate all variogram-related inputs
validate_variogram_parameters(variogram_parameters, fit_parameters, optimization_settings)

best_fit_variogram, initial_fit, optimized_fit = optimize_variogram(
    lag_counts, lags, gamma_h, variogram_parameters, optimization_settings
)
best_fit_variogram

survey = Survey(init_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/initialization_config.yml" ,
                survey_year_config_path="C:/Users/Brandyn/Documents/GitHub/echopop/config_files/survey_year_2019_config.yml")
survey.transect_analysis()
survey.fit_variogram(jacobian_approx="central")
self = survey
settings_dict = {"variable": "biomass_density",
                 "stratum_name": "stratum_num",
                 "kriging_parameters": {
                    "longitude_reference": -124.78338,
                    "latitude_offset": 45.0,
                    "longitude_offset": -124.78338,
                 }}

# Prepare the transect data
# ---- Create a copy of the transect dictionary
transect_input = copy.deepcopy(self.analysis["transect"])
# ---- Edit the transect data
transect_data = edit_transect_columns(transect_input, settings_dict)
isobath_df = self.input["statistics"]["kriging"]["isobath_200m_df"]
transect_data, _, _ = transform_geometry(transect_data, isobath_df, settings_dict)

# Create a copy of the existing variogram settings
variogram_parameters = self.input["statistics"]["variogram"]["model_config"].copy()
n_lags = 30
lag_resolution = variogram_parameters["lag_resolution"]
# ---- Number of lags
variogram_parameters["n_lags"] = 30
# ---- Azimuth range
variogram_parameters["azimuth_range"] = 360
# ---- Force lag-0
variogram_parameters["force_lag_zero"] = True
# Compute the lag distances
distance_lags = np.arange(1, n_lags) * lag_resolution
# ---- Add to the `variogram_parameters` dictionary
variogram_parameters["distance_lags"] = distance_lags
# ---- Update the max range parameter, if necessary
max_range = lag_resolution * n_lags
lags, gamma_h, lag_counts, lag_covariance = empirical_variogram(
    transect_data, variogram_parameters, {"variable": "biomass_density"}
)
fig1, ax1 = plt.subplots()
ax1.imshow(np.random.random((5,5)))
plt.show()
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
    return fig

fig = plot_variogram()

plot_variogram()
plt.show()

def plot_variogram():
    s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
    plt.figure(figsize=(6,4))
    plt.scatter(lags, gamma_h, s=s(lag_counts))
    plt.xlabel('Lag Distance')
    plt.ylabel('Semivariance')
    plt.title('Empirical Variogram')
    plt.grid(True)
    plt.tight_layout()

    # Embed plot into tkinter GUI
    canvas = FigureCanvasTkAgg(plt.gcf(), master=plot_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Optionally add toolbar for plot navigation
    toolbar = NavigationToolbar2Tk(canvas, plot_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

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


param_key = {
    "nugget": dict(name="Nugget", widget="entry"),
    "sill": dict(name="Sill", widget="entry"),
    "correlation_range": dict(name="Correlation range", widget="entry"),
    "decay_power": dict(name="Decay power", widget="entry"),
    "hole_effect_range": dict(name="Hole effect range", widget="entry"),
    "enhance_semivariance": dict(name="Enhance semivariance", widget="checkbox")
}

optimization_key = {
    "max_fun_evaluations": dict(name="Maximum function evaluations", widget="entry"),
    "cost_fun_tolerance": dict(name="Cost function tolerance", widget="entry"),
    "solution_tolerance": dict(name="Solution tolerance", widget="entry"),
    "gradient_tolerance": dict(name="Gradient tolerance", widget="entry"),
    "finite_step_size": dict(name="Finite differences step size", widget="entry"),
    "trust_region_solver": dict(name="Trust Region solver method", widget="dropdown"),
    "x_scale": dict(name="Characteristic x-scale", widget=["checkbox", "entry"]),
    "jacobian_approx": dict(name="Jacobian approximation method", widget="dropdown"),
}

optimization_parameters = {
    "max_fun_evaluations": max_fun_evaluations,
    "cost_fun_tolerance": cost_fun_tolerance,
    "solution_tolerance": solution_tolerance,
    "gradient_tolerance": gradient_tolerance,
    "finite_step_size": finite_step_size,
    "trust_region_solver": trust_region_solver,
    "x_scale": x_scale,
    "jacobian_approx": jacobian_approx,
}

def get_optimization_defaults(optimization_parameters):
    # ----


def on_selection(event):
    """
    """
    selected_model = combobox.get()
    get_variogram_model(selected_model)

def get_variogram_model(selected_model = "Bessel-exponential"):
    """
    """

    # Get the corresponding value from the `WIDGET_MAP`
    model = VARIOGRAM_MODEL_WIDGET_MAP.get(selected_model)

    # Begin processing
    if isinstance(model, list):
        print(f"Processing composite model: {selected_model}")
    else:
        print(f"Processing model: {selected_model}")

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
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Create the main application window
root = tk.Tk()
root.title("Fit Variogram GUI")

# Create a main frame to hold all components
main_frame = ttk.Frame(root)
main_frame.pack(fill=tk.BOTH, expand=True)

# Create a frame for the left side (tabs)
left_frame = ttk.Frame(main_frame)
left_frame.grid(row=0, column=0, sticky=tk.NSEW, padx=10, pady=10)

# Create a frame for the right side (plot panel)
right_frame = ttk.Frame(main_frame)
right_frame.grid(row=0, column=1, sticky=tk.NSEW, padx=10, pady=10)

# Create tabs using ttk.Notebook in the left frame
notebook = ttk.Notebook(left_frame)
notebook.pack(fill=tk.BOTH, expand=True)

# Create frames for each tab
variogram_frame = ttk.Frame(notebook)
optimization_frame = ttk.Frame(notebook)

# Add tabs to the notebook
notebook.add(variogram_frame, text='Variogram Parameters')
notebook.add(optimization_frame, text='Optimization Parameters')

# Dropdown widget
combobox = ttk.Combobox(left_frame, values=list(VARIOGRAM_MODEL_WIDGET_MAP.keys()), state="readonly", width=20)
combobox.set("Bessel-exponential")
combobox.pack(padx=10, pady=5)





# Function to handle x-scale selection
def toggle_x_scale():
    if x_scale_var.get():
        x_scale_entry.config(state=tk.DISABLED)
        x_scale_entry.delete(0, tk.END)  # Clear the entry field when disabled
    else:
        x_scale_entry.config(state=tk.NORMAL)
        


# Function to update fit parameters based on selected model
def update_fit_parameters(event):

    param_key = {
        "nugget": "Nugget",
        "sill": "Sill",
        "correlation_range": "Correlation range",
        "decay_power": "Decay power",
        "hole_effect_range": "Hole effect range"
    }
    selected_model = combobox.get()
    print(f"Selected Model: {selected_model}")

    # Clear existing widgets in the variogram frame
    for widget in variogram_frame.winfo_children():
        widget.destroy()

    # Generate checkboxes and entry fields based on selected model
    if selected_model in VARIOGRAM_MODEL_WIDGET_MAP:
        fit_params = VARIOGRAM_MODEL_WIDGET_MAP[selected_model]
        function_arguments, _ = get_variogram_arguments(fit_params)

        for param in function_arguments.values():
            if param.name != "distance_lags":
                arg = param.name
                # ---- get the name from the key
                arg_name = param_key[arg]
                var = tk.BooleanVar(value=True)
                checkbox = ttk.Checkbutton(variogram_frame, text=arg_name, variable=var)
                checkbox.pack(padx=5, pady=5)

                # Add an entry field next to each checkbox
                entry = ttk.Entry(variogram_frame, width=10)
                entry.pack(padx=5, pady=5)
    else:
        print(f"No fit parameters defined for model: {selected_model}")

# Bind selection event to update fit parameters
combobox.bind("<<ComboboxSelected>>", update_fit_parameters)

# Function to create the optimization parameters frame
def create_optimization_frame():
    for field, var_name in optimization_fields:
        label = ttk.Label(optimization_frame, text=field)
        label.pack(padx=10, pady=2)

        if var_name == "trust_region_solver":
            combobox = ttk.Combobox(optimization_frame, values=["Exact", "Base"], state="readonly")
            combobox.pack(padx=10, pady=2)
        elif var_name == "jacobian_approx":
            combobox = ttk.Combobox(optimization_frame, values=["Forward", "Central"], state="readonly")
            combobox.pack(padx=10, pady=2)
        else:
            entry = ttk.Entry(optimization_frame)
            entry.pack(padx=10, pady=2)

# Initialize optimization frame
# create_optimization_frame()
# Create widgets for each optimization parameter
for field, var_name in optimization_fields:
    label = ttk.Label(optimization_frame, text=field)
    label.pack(padx=10, pady=2)

    if var_name == "x_scale":
        x_scale_var = tk.BooleanVar(value=True)  # Default checked
        x_scale_checkbox = ttk.Checkbutton(optimization_frame, text="Use Jacobian", variable=x_scale_var,
                                           command=toggle_x_scale)
        x_scale_checkbox.pack(padx=10, pady=2)

        x_scale_entry = ttk.Entry(optimization_frame, state=tk.DISABLED)
        x_scale_entry.pack(padx=10, pady=2)

        # Initial setup based on default value
        toggle_x_scale()  # Call toggle function initially

    else:
        entry = ttk.Entry(optimization_frame)
        entry.pack(padx=10, pady=2)
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

def plot_variogram():
    s = lambda x : (((x-x.min())/float(x.max()-x.min())+1)*10)**2
    x = lags; y = gamma_h
    # plt.figure(figsize=(8,6))
    # plt.scatter(x, y, s=s(lag_counts))
    # plt.xlabel('Lag Distance')
    # plt.ylabel('Semivariance')
    # plt.title('Empirical Variogram')
    # plt.grid(True)
    # plt.tight_layout()

    fig, ax = plt.subplots(figsize=(8, 6))
    scatter = ax.scatter(x, y, s=s(lag_counts))
    ax.set_xlabel('Lag Distance')
    ax.set_ylabel('Semivariance')
    ax.set_title('Empirical Variogram')
    ax.grid(True)
    plt.tight_layout()

    # Embed plot into tkinter GUI
    canvas = FigureCanvasTkAgg(plt.gcf(), master=right_frame)
    canvas.draw()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

    # Optionally add toolbar for plot navigation
    toolbar = NavigationToolbar2Tk(canvas, right_frame)
    toolbar.update()
    canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

     # Create a list of annotated points
    annot = ax.annotate("", xy=(0,0), xytext=(10,10), textcoords="offset points",
                            bbox=dict(boxstyle="round", fc="white", ec="black", alpha=1.0),
                            arrowprops=dict(arrowstyle="->", color="black", lw=1, linestyle="-"))
    annot.set_visible(False)

    def update_annot(ind):
        pos = scatter.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = f"Lag={pos[0]:.2f}\nγ={pos[1]:.2f}\nCount={lag_counts[ind['ind'][0]]}"
        annot.set_text(text)
        annot.get_bbox_patch().set_alpha(1.0)

    def on_hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = scatter.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                fig.canvas.draw_idle()
            else:
                if vis:
                    annot.set_visible(False)
                    fig.canvas.draw_idle()

    fig.canvas.mpl_connect("motion_notify_event", on_hover)

# Button to compute empirical variogram
button_compute = ttk.Button(left_frame, text="Compute Empirical Variogram", command=plot_variogram)
button_compute.pack(padx=10, pady=10)

# Start the main event loop
root.mainloop()