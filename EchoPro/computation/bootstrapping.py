from .kriging import krig_type_dict
from ..data_loader import KrigingMesh
from .kriging import Kriging
from typing import Optional, Tuple, List
import math
import numpy as np
import pandas as pd


class Bootstrapping:
    """
    A Class that performs bootstrapping for data that has been
    Kriged and data that has not been Kriged.

    Parameters
    ----------
    survey : Survey
        An initialized Survey object.
    """

    def __init__(self, survey):

        # initialize class survey object
        self.survey = survey

    def _check_kriging_params(self, kriging_params: dict) -> Tuple[dict, KrigingMesh]:
        """
        Ensure that the input ``kriging_params`` contains all
        necessary parameters and they are the appropriate type.
        Additionally, assigns the ``krig_mesh`` parameter to
        its own variable and removes it from ``kriging_params``.

        Parameters
        ----------
        kriging_params: dict
            All parameters needed for Kriging in addition to the
            ``krig_mesh`` parameter.

        Returns
        -------
        kriging_params: dict
            Dictionary with all necessary Kriging parameters and
            ``krig_mesh`` removed
        krig_mesh_obj: KrigingMesh
            A KrigingMesh object obtained from the ``krig_mesh``
            key in ``kriging_params``
        """

        # make sure an empty dict was not given
        if kriging_params is None:
            raise ValueError("If Kriging should be run, the input ``kriging_params`` must be provided!")

        # ensure that krig_mesh was provided and is of the appropriate type
        if ("krig_mesh" not in kriging_params.keys()) or (not isinstance(kriging_params["krig_mesh"], KrigingMesh)):
            raise TypeError("The parameter ``krig_mesh`` was not provided or it is not of type KrigingMesh!")

        # get krig_mesh object and remove it from kriging_params
        krig_mesh_obj = kriging_params["krig_mesh"]
        del kriging_params["krig_mesh"]

        # make sure all kriging parameters are included
        if set(krig_type_dict.keys()).difference(set(kriging_params.keys())):
            raise ValueError("Some required Kriging parameters where not provided!")

        # check that all types are correct for params
        for key, val in kriging_params.items():
            expected_type = krig_type_dict.get(key)
            if not isinstance(val, expected_type):
                raise TypeError(f"The Kriging parameter {key} is not of type {expected_type}")

        # make sure that the krig_mesh object refers to the appropriate survey object
        if id(krig_mesh_obj.survey) != id(self.survey):
            raise RuntimeError("The krig_mesh Survey object and the Survey object used to "
                               "initialize bootstrapping must be the same!")

        return kriging_params, krig_mesh_obj

    def _get_results_for_no_kriging(self) -> List[float]:
        """
        Obtains the total areal biomass density and associated
        mean Jolly-Hampton CV value for a bootstrapping iteration.

        Returns
        -------
        tot_bio_mass_no_kriging: float
            The total areal biomass density for data that has not
            been Kriged
        CV_JH_mean_no_kriging: float
            The mean Jolly-Hampton CV value for the data that has not
            been Kriged
        """

        # calculate total biomass density
        tot_bio_mass_no_kriging = self.survey.bio_calc.final_biomass_table["areal_biomass_density_adult"].sum()

        # perform CV analysis on data
        CV_JH_mean_no_kriging = self.survey.run_cv_analysis(kriged_data=False)

        return [tot_bio_mass_no_kriging, CV_JH_mean_no_kriging]

    def _get_results_for_kriging(self, krig_mesh_obj: KrigingMesh,
                                 krig: Kriging) -> List[float]:
        """
        Obtains the total Kriged biomass estimate and associated
        mean Jolly-Hampton CV value for a bootstrapping iteration.

        Parameters
        ----------
        krig_mesh_obj: KrigingMesh
            An initialized ``KrigingMesh`` object used to run Kriging
        krig: Kriging
            An initialized ``Kriging``object

        Returns
        -------
        tot_bio_mass_kriging: float
            The total biomass estimate for data that has
            been Kriged
        CV_JH_mean_kriging: float
            The mean Jolly-Hampton CV value for the data that has
            been Kriged
        """

        # TODO: we should make it so that the user can provide their
        #  own transformation for the below meshes

        # apply transect mesh transformation
        krig_mesh_obj.apply_coordinate_transformation(coord_type='transect')

        # apply full mesh transformation
        # TODO: they way this is setup is inefficient because we need
        #  to transform the full mesh every time, can we avoid this?
        krig_mesh_obj.apply_coordinate_transformation(coord_type='mesh')

        # get Kriged biomass estimate
        krig.run_biomass_kriging(krig_mesh_obj)

        # calculate the total Kriged biomass density
        tot_bio_mass_kriging = self.survey.bio_calc.krig_results_gdf.biomass.sum()

        # perform CV analysis on Kriged data
        CV_JH_mean_kriging = self.survey.run_cv_analysis(kriged_data=True)

        return [tot_bio_mass_kriging, CV_JH_mean_kriging]

    def run_bootstrapping(self, run_kriging: bool = False,
                          kriging_params: Optional[dict] = None,
                          removal_percentage: float = 50.0,
                          num_iterations: int = 10,
                          seed: Optional[int] = None) -> pd.DataFrame:
        """
        A routine for running bootstrapping on a reduced number of
        transects within a survey. The bootstrapping performed can
        be completed for data that has been Kriged and data that
        has no Kriging. Each bootstrapping iteration is obtained
        by randomly selecting a subset of transects based on the
        ``removal_percentage`` input.

        Parameters
        ----------
        run_kriging: bool
            If True, Kriging will be run (this requires that one provide the
            input ``kriging_params``), otherwise Kriging will not be run
        kriging_params: dict or None
            All parameters needed to initialize the kriging routine via
            ``survey.get_kriging()`` in addition to the parameter ``krig_mesh``,
            which is an initialized ``KrigingMesh`` object
        removal_percentage: float
            The percentage of transects that should be removed.
        num_iterations: int
            The number of bootstrapping iterations to run
        seed: int or None
            The seed for the random number generator used to select the transects. If
            no seed is provided the random number generator will not be seeded.

        Returns
        -------
        boot_final_results: pd.Dataframe
            A Dataframe corresponding to the results for each bootstrap iteration.

        Notes
        -----
        If one wants to remove 30% of the transects and then run bootstrapping,
        ``removal_percentage`` should be set to ``30.0``.

        If ``removal_percentage`` multiplied by the number of unique transects in
        ``survey.nasc_df`` is not an integer, then the ceiling of the result will
        be used for the number of transects to be removed.

        If ``run_kriging=False`` then the Dataframe will have the columns
        ``tot_biomass_no_kriging`` and ``CV_JH_mean_no_kriging``, which correspond to
        the total biomass density with no Kriging and the mean of the Jolly-Hampton CV value
        for the associated data, respectively.

        If ``run_kriging=True`` then the Dataframe will have the same columns as when
        ``run_kriging=False`` in addition to the columns ``tot_biomass_kriging`` and
        ``CV_JH_mean_kriging``, which correspond to the total biomass density produced by the Kriged
        data and the mean of the Jolly-Hampton CV value for the associated data, respectively.

        The number of Jolly-Hampton realizations is determined by an initialization parameter,
        see ``Survey.run_cv_analysis`` for more details.
        """

        if run_kriging:

            # check kriging_params input and obtain KrigingMesh object
            kriging_params, krig_mesh_obj = self._check_kriging_params(kriging_params)

            # initalize kriging routine
            krig = self.survey.get_kriging(kriging_params)

        # get all unique transects in nasc_df
        unique_transects = self.survey.nasc_df.index.unique().values

        # determine the number of transects that should be randomly selected
        num_sel_transects = math.floor(len(unique_transects) * (1.0 - removal_percentage / 100.0))

        # initialize the random number generator object
        rng = np.random.default_rng(seed)

        vals_to_keep = []
        for iteration in range(num_iterations):

            # randomly select transects without replacement
            selected_transects = list(rng.choice(unique_transects, num_sel_transects, replace=False))

            # compute the areal biomass density for the subset of transects
            self.survey.compute_biomass_density(selected_transects=selected_transects)

            # collect total biomass and associated JH CV value for data without Kriging
            vals_to_keep.append(self._get_results_for_no_kriging())

            if run_kriging:
                # collect total biomass and associated JH CV value for data with Kriging
                vals_to_keep[iteration] += self._get_results_for_kriging(krig_mesh_obj, krig)

        # set DataFrame column names
        col_names = ["tot_biomass_no_kriging", "CV_JH_mean_no_kriging"]
        if run_kriging:
            col_names += ["tot_biomass_kriging", "CV_JH_mean_kriging"]

        # construct final DataFrame holding bootstrapping results
        boot_final_results = pd.DataFrame(vals_to_keep,
                                          columns=col_names, index=range(1, num_iterations+1),
                                          dtype=np.float64)
        boot_final_results.index.name = "iteration"

        return boot_final_results
