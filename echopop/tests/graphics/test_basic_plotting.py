import matplotlib
import pytest

from echopop import graphics as egra

# Change to non-interactive GUI
matplotlib.use("Agg")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_kriged_mesh_types(plotting_mesh_data):

    # Hexbin (implicit)
    egra.plot_kriged_mesh(plotting_mesh_data, "biomass")

    # Hexbin (explicit)
    egra.plot_kriged_mesh(plotting_mesh_data, "biomass", plot_type="hexbin")

    # Scatter
    egra.plot_kriged_mesh(plotting_mesh_data, "biomass", plot_type="scatter")

    # Pseudocolor mesh
    egra.plot_kriged_mesh(plotting_mesh_data, "biomass", plot_type="pcolormesh")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_kriged_mesh_invalid(plotting_mesh_data):
    with pytest.raises(ValueError):
        egra.plot_kriged_mesh(plotting_mesh_data, "biomass", plot_type="not_a_type")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_transect_map(plotting_transect_data):
    egra.plot_transect_map(plotting_transect_data, "biomass")


@pytest.mark.filterwarnings("ignore::UserWarning")
def test_plot_age_length_heatmap(plotting_heatmap_data):
    egra.plot_age_length_heatmap(plotting_heatmap_data)
