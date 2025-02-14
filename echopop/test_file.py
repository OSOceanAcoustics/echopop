from echopop.survey import Survey
from echopop.extensions import inversion

# Initialize the `Survey` object with the inversion parameterization and dataset
inversion_configuration_filepath = "C:/Users/Brandyn/Documents/GitHub/echopop/config_files/inversion_processing.yml"
survey = Survey.initialize_inversion(inversion_configuration_filepath)

# Run inversion algorithm (see `inversion_processing.yml` for parameterization)
survey.invert_population(subset_dataset={"transect_min": 1, "transect_max": 60})
# ---- Parse results
survey.results["inversion"]

biomass_df = survey.results["inversion"]['transect_df']

import matplotlib.pyplot as plt
import numpy as np
plt.hist(np.log10(biomass_df["biomass_density"] + 1))
plt.show()

# Plot the transect data
survey.plot(kind="transect", variable="biomass_density")

# Optimize normalized variogram to the inverted population estimates
survey.fit_inversion_variogram(
    model="exponential",
    n_lags=60,
    initialize_variogram = {
        "nugget": dict(min=0.0, max=0.7, value=0.1, vary=True),
        "sill": dict(min=0.5, max=1.5, value=0.9, vary=True),
        "correlation_range": dict(min=1e-3, max=1.0, value=0.1, vary=True)
    },
    optimization_parameters={
        "max_fun_evaluations": 5000,
        "finite_step_size": 1e-10,
        "jacobian_approx": "central",
    }
)
# ---- Parse results
survey.results["variogram"]

# Krige the inverted population estimates using the best-fit variogram
survey.inversion_kriging_analysis()
# ---- Pares results
survey.results["kriging"]

# Plot the kriged estimates
survey.plot(kind="mesh", variable="biomass", plot_type="pcolormesh")

