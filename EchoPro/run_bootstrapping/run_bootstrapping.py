from ..global_vars import krig_type_dict
from ..kriging_mesh import KrigingMesh
from typing import Optional
import math
import numpy as np
import pandas as pd


class Bootstrapping:
    """
    A Class that performs bootstrapping. This involves processing
    acoustic and biological data, computation of CV analysis,
    and Kriging.  TODO: modify this doc

    Parameters
    ----------
    survey : Survey
        An initialized Survey object. Note that any change to
        ``self.survey`` will also change this object.
    """

    def __init__(self, survey):

        self.survey = survey

    def run_bootstrapping(self, run_kriging: bool = False,
                          kriging_params: Optional[dict] = None,
                          removal_percentage: float = 50.0,
                          num_iterations: int = 10,
                          seed: Optional[int] = None):
        """



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

        Notes
        -----
        If one wants to remove 30% of the transects and then run bootstrapping,
        then ``removal_percentage`` should be set to ``30.0``.

        If ``removal_percentage`` multiplied by the number of unique transects in
        ``self.nasc_df`` is it not an integer, then the ceiling of the result will
        be used for the number of transects to be removed.
        """

        if run_kriging:

            # TODO: move all of this to a function

            if kriging_params is None:
                raise ValueError("If Kriging should be run, the input ``kriging_params`` must be provided!")

            if ("krig_mesh" not in kriging_params.keys()) or (not isinstance(kriging_params["krig_mesh"], KrigingMesh)):
                raise TypeError("The parameter ``krig_mesh`` was not provided or it is not of type KrigingMesh!")

            # get krig_mesh object and remove it from kriging_params
            krig_mesh_obj = kriging_params["krig_mesh"]
            del kriging_params["krig_mesh"]

            # make sure all kriging parameters are included
            if set(krig_type_dict.keys()).difference(set(kriging_params.keys())):
                raise ValueError("Some required parameters where not provided!")

            # check that all types are correct for params
            for key, val in kriging_params.items():
                expected_type = krig_type_dict.get(key)
                if not isinstance(val, expected_type):
                    raise TypeError(f"The Kriging parameter {key} is not of type {expected_type}")

            # TODO: make sure that kriging_params["krig_mesh"].survey and
            #  `self.survey` are the same objects, error out if they are not

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

            self.survey.compute_biomass_density(selected_transects=selected_transects)

            tot_bio_mass_no_kriging = self.survey.bio_calc.final_biomass_table["normalized_biomass_density"].sum()

            CV_JH_mean_no_kriging = self.survey.run_cv_analysis(kriged_data=False)

            if run_kriging:

                # TODO: we should make it so that the user can provide their
                #  own transformation for the below meshes

                # apply transect mesh transformation
                krig_mesh_obj.apply_coordinate_transformation(coord_type='transect')

                # apply full mesh transformation
                # TODO: they way this is setup is inefficient because we need
                #  to transform the full mesh every time, can we avoid this?
                krig_mesh_obj.apply_coordinate_transformation(coord_type='mesh')

                krig.run_biomass_kriging(krig_mesh_obj)

                tot_bio_mass_kriging = self.survey.bio_calc.krig_results_gdf.krig_biomass_vals.sum()
                CV_JH_mean_kriging = self.survey.run_cv_analysis(kriged_data=True)

                vals_to_keep.append([tot_bio_mass_no_kriging, CV_JH_mean_no_kriging,
                                     tot_bio_mass_kriging, CV_JH_mean_kriging])

            else:

                vals_to_keep.append([tot_bio_mass_no_kriging, CV_JH_mean_no_kriging])

        if run_kriging:
            col_names = ["tot_biomass_no_kriging", "CV_JH_mean_no_kriging",
                         "tot_biomass_kriging", "CV_JH_mean_kriging"]
        else:
            col_names = ["tot_biomass_no_kriging", "CV_JH_mean_no_kriging"]

        boot_final_results = pd.DataFrame(vals_to_keep,
                                          columns=col_names, index=range(1, num_iterations+1),
                                          dtype=np.float64)

        boot_final_results.index.name = "iteration"

        return boot_final_results
