import ipywidgets
import numpy as np
from scipy import special
import ipywidgets as widgets
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import inspect
from ..numba_modules import nb_subtract_outer, nb_dis_vec, nb_diff_sqrd
from typing import Dict, Tuple, Callable

# default bounds for fitting the semi-variogram model
bnds_default = {'LB': {'sill': 0.0, 'ls': 0.0,
                       'exp_pow': 0.0, 'ls_hole_eff': 0.0, 'nugget': 0.0},
                'UB': {'sill': np.inf, 'ls': np.inf,
                       'exp_pow': np.inf, 'ls_hole_eff': 1e-13, 'nugget': 1e-13}}


class SemiVariogram:
    """
    This class contains a routine that calculates a
    normalized semi-variogram and routines for
    obtaining the best semi-variogram model for the
    estimated semi-variogram.

    Parameters
    ----------
    x : np.ndarray
        A 1D array representing the x coordinate for
        the semi-variogram calculation
    y : np.ndarray
        A 1D array representing the y coordinate for
        the semi-variogram calculation
    field : np.ndarray
        A 1D array representing the field for the
        semi-variogram calculation (e.g. biomass density)
    lag_res : float
        The spacing between lag centers
    nlag : int
        The total number of lag centers

    Notes
    -----
    The variables ``lag_res`` and ``nlag`` determine
    the center of bins used to calculate the
    semi-variogram. For example, if `lag_res=3`` and
    ``nlag=0.1`` then we would obtain the bin centers
    ``[0.0, 0.1, 0.2]``.
    """

    def __init__(self, x: np.ndarray, y: np.ndarray,
                 field: np.ndarray, lag_res: float, nlag: int):

        # input data
        self.x = x
        self.y = y
        self.field = field
        self.nlag = nlag
        self.lag_res = lag_res

        # bins provided to calculate_semi_variogram
        self._center_bins = lag_res * np.arange(nlag)

        # normalized semi-variogram values
        self.gamma_normalized = None

        # widget variables
        self._model_drop = None
        self._sill_box = None
        self._len_scale_box = None
        self._exp_pow_box = None
        self._ls_hole_eff_box = None
        self._nugget_box = None
        self._apply_mod_toggle = None
        self._lsq_toggle = None
        self._full_params = None

    def calculate_semi_variogram(self, center_bins: np.ndarray = None) -> None:
        """
        Calculates the semi-variogram normalized by the
        standard deviation of the head multiplied by
        the standard deviation of the tail for each
        lag. This calculation assumes that the mesh
        points (x,y) are isotropic and the search
        area is omnidirectional.

        Parameters
        ----------
        center_bins: np.ndarray
            A 1D array representing the center of the bins used
            to calculate the semi-variogram

        Notes
        -----
        The calculated normalized semi-variogram values are
        stored in the class variable ``gamma_normalized``.
        """

        # check center_bins and assign it, if necessary
        if isinstance(center_bins, np.ndarray):
            if center_bins.ndim == 1:
                center_bins = np.sort(center_bins)
            else:
                ValueError("center_bins must be 1-dimensional!")
        elif center_bins is None:
            center_bins = self._center_bins
        else:
            raise ValueError("center_bins must be an np.ndarray or None!")

        # get upper triangular indices
        i_up_ind, j_up_ind = np.triu_indices(len(self.x), k=1)

        # get the upper triangular mesh of differences for x/y
        x_diff = nb_subtract_outer(self.x, self.x)[i_up_ind, j_up_ind]
        y_diff = nb_subtract_outer(self.y, self.y)[i_up_ind, j_up_ind]

        # find the distance between points
        dis = nb_dis_vec(x_diff, y_diff)

        # construct head and tail field values
        field_rep = np.tile(self.field, (len(self.field), 1))
        field_head = field_rep[j_up_ind, i_up_ind]
        field_tail = field_rep[i_up_ind, j_up_ind]

        # get the squared difference of the head and the tail
        field_diff_sqrd = nb_diff_sqrd(field_head, field_tail)

        self.gamma_normalized = np.empty(len(center_bins), dtype=np.float64)
        for i in range(len(center_bins)):
            # get indices of distances that are in the lag
            if i == 0:
                bin_add_m = (center_bins[i + 1] - center_bins[i]) / 2.0
                bin_add_p = (center_bins[i + 1] - center_bins[i]) / 2.0
            elif i == len(center_bins) - 1:
                bin_add_m = (center_bins[i] - center_bins[i - 1]) / 2.0
                bin_add_p = (center_bins[i] - center_bins[i - 1]) / 2.0
            else:
                bin_add_p = (center_bins[i + 1] - center_bins[i]) / 2.0
                bin_add_m = (center_bins[i] - center_bins[i - 1]) / 2.0

            ind_in_lag = np.argwhere((center_bins[i] - bin_add_m <= dis)
                                     & (dis < center_bins[i] + bin_add_p))

            # calculate the semi-variogram value
            gamma = 0.5 * np.mean(field_diff_sqrd[ind_in_lag])

            # normalize gamma by the standard deviation of the head
            # multiplied by the standard deviation of the tail
            std_head = np.std(field_head[ind_in_lag])
            std_tail = np.std(field_tail[ind_in_lag])
            self.gamma_normalized[i] = gamma / (std_head * std_tail)

    def _create_widgets(self):
        """
        Constructs all widgets needed to create the
        ``view_semi_variogram`` interactive output.
        Additionally, stores all of these widgets
        as class variables.
        """

        # widget to define semi-variogram model
        self._model_drop = widgets.Dropdown(options=['Exponential', 'Generalized exponential-Bessel'],
                                            value='Generalized exponential-Bessel',
                                            description='Semi-variogram model',
                                            disabled=False, style={'description_width': 'initial'})

        # widgets that accept/display semi-variogram model parameters
        self._sill_box = widgets.FloatText(value=1.0,
                                           description='Sill', disabled=False)
        self._len_scale_box = widgets.FloatText(value=1.0,
                                                description='Length Scale',
                                                disabled=False,
                                                style={'description_width': 'initial'})
        self._exp_pow_box = widgets.FloatText(value=1.0,
                                              description='Exponential power',
                                              disabled=False,
                                              style={'description_width': 'initial'})
        self._ls_hole_eff_box = widgets.FloatText(value=1.0,
                                                  description='Length scale hole effect',
                                                  disabled=False,
                                                  style={'description_width': 'initial'})
        self._nugget_box = widgets.FloatText(value=0.0,
                                             description='Nugget', disabled=False)

        # collect all widgets that define the semi-variogram model parameters
        self._full_params = {'sill': self._sill_box, 'ls': self._len_scale_box,
                             'exp_pow': self._exp_pow_box, 'ls_hole_eff': self._ls_hole_eff_box,
                             'nugget': self._nugget_box}

        # widget that triggers the Least Squares fit
        self._lsq_toggle = widgets.ToggleButton(value=False,
                                                description='Run Least Squares fit',
                                                disabled=False, button_style='info',
                                                tooltip='Description',
                                                style={'description_width': 'initial'})

        # widget to apply the model
        self._apply_mod_toggle = widgets.ToggleButton(value=False,
                                                      description='Apply model',
                                                      disabled=False, button_style='info',
                                                      tooltip='Description',
                                                      icon='line-chart',
                                                      style={'description_width': 'initial'})

    @staticmethod
    def generalized_exp_bessel(lag_vec: np.ndarray, sill: float, ls: float,
                               exp_pow: float, ls_hole_eff: float,
                               nugget: float) -> np.ndarray:
        """
        Applies a generalized exponential bessel function
        to the provided input.

        Parameters
        ----------
        lag_vec : np.ndarray
            A 1D array where the elements correspond to the lag
        sill : float
            Sill of model
        ls : float
            Length scale for the main lobe
        exp_pow : float
            Power for the exponential arguments
        ls_hole_eff : float
            Length scale for hole effect
        nugget : float
            Nugget effect

        Returns
        -------
        np.ndarray
            1D array of values obtained from evaluating the function
        """

        return (sill - nugget)*(1.0 - np.exp(-(lag_vec/ls)**exp_pow)*special.j0(ls_hole_eff*lag_vec)) + nugget

    @staticmethod
    def exponential(lag_vec: np.ndarray, sill: float, ls: float,
                    nugget: float) -> np.ndarray:
        """
        Applies an function to the provided input.

        Parameters
        ----------
        lag_vec : np.ndarray
            A 1D array where the elements correspond to the lag
        sill : float
            Sill of model
        ls : float
            Length scale for the main lobe
        nugget : float
            Nugget effect

        Returns
        -------
        np.ndarray
            1D array of values obtained from evaluating the function
        """

        return (sill - nugget)*(1.0 - np.exp(-(lag_vec/ls))) + nugget

    def get_model_info_from_str(self, model: str) -> Tuple[Callable, tuple]:
        """
        Selects the appropriate model and model
        kwargs needed, based on the provided model
        name.

        Parameters
        ----------
        model : str
            The name of the model

        Returns
        -------
        Tuple[Callable, tuple]
            The first element is a callable representing
            the model and the second is the required kwargs
            for that model
        """

        # select corresponding function
        if model == 'Generalized exponential-Bessel':
            func = self.generalized_exp_bessel
        elif model == 'Exponential':
            func = self.exponential

        # get required kwargs for func
        model_keys = tuple(inspect.signature(func).parameters)

        return func, model_keys

    def _evaluate_model(self, params: dict, model: str) -> np.ndarray:
        """
        Selects the appropriate model function and then
        evaluates it using the provided parameters.

        Parameters
        ----------
        params : dict
            A dictionary of all parameters needed to
            evaluate the model
        model : str
            The name of the model to evaluate

        Returns
        -------
        np.ndarray
            1D array of values obtained from evaluating the model
        """

        func, model_keys = self.get_model_info_from_str(model)

        # select only those keys in params that are needed by func
        sub_dict = dict((k, params[k]) for k in model_keys if k in params)

        # evaluate model
        return func(*tuple(sub_dict.values()))

    def apply_lsq_fit(self, model: str) -> Dict[str, float]:
        """
        Obtains the set of model parameters that minimize the sum
        of the squared residuals. The residual is defined as the
        difference between the normalized semi-variogram and
        the semi-variogram model selected.

        Parameters
        ----------
        model : str
            The name of the model to use

        Returns
        -------
        Dict[str, float]
            Dictionary with the model parameters as keys and
            the parameters that minimize the residual as values

        Notes
        -----
        This function uses a pre-defined set of boundaries for
        the fitting.
        # TODO: allow for user specification of boundaries
        """

        func, model_keys = self.get_model_info_from_str(model)

        # create lower and upper bounds for the fitting parameters
        low_bnds = [bnds_default['LB'][key] for key in model_keys if key in bnds_default['LB']]
        up_bnds = [bnds_default['UB'][key] for key in model_keys if key in bnds_default['UB']]

        # fit the model to the normalized semi-variogram values
        fit_params, _ = curve_fit(func, self._center_bins, self.gamma_normalized,
                                  bounds=(low_bnds, up_bnds))

        # TODO: Open up an issue for the differences in lsq against Matlab version

        # collect model parameters
        return dict(zip(model_keys[1:], fit_params))

    def _plot_model(self, all_inputs: dict, model: str,
                   ax: plt.axes) -> plt.axes:
        """
        Plots a semi-variogram model using the
        provided parameters on the axis ``ax``.

        Parameters
        ----------
        all_inputs : dict
            A dictionary of all parameters needed to
            evaluate the model
        model : str
            The name of the model to evaluate
        ax : plt.axes
            The axis to plot on

        Returns
        -------
        plt.axes
            Input axis with the model plotted on it
        """

        # create lag values to plot along
        all_inputs['lag_vec'] = np.linspace(np.min(self._center_bins),
                                            np.max(self._center_bins), 100)

        model_vals = self._evaluate_model(all_inputs, model)

        ax.plot(all_inputs['lag_vec'], model_vals, 'r-')

        return ax

    def _manage_plots(self, model: str, lsq_fit: bool, apply_model: bool) -> None:
        """
        Manages the construction of the figure to plot on,
        plotting of the normalized semi-variogram values,
        Least Squares fitting, and plotting of the model.

        Parameters
        ----------
        model : str
            The name of the model to evaluate/plot
        lsq_fit : bool
            If True, Least Squares fitting is applied,
            otherwise nothing is done
        apply_model : bool
            If True, plots the model for a given set of parameters,
            otherwise does nothing
        """

        plt.close(999)  # close the figure, if it exists

        # create a figure to hold the plotted data
        fig = plt.figure(999, figsize=(8, 4), dpi=100)
        fig.canvas.manager.set_window_title('')
        plt.close(999)

        ax1 = fig.add_subplot()  # create axis to plot on

        # plot the normalized semi-variogram values
        ax1.plot(self._center_bins, self.gamma_normalized, 'bo')
        ax1.set_ylabel("gamma normalized")
        ax1.set_xlabel("lag")

        # apply Least Squares, plot the data, and update the parameters
        if lsq_fit:
            lsq_fit_vals = self.apply_lsq_fit(model)
            self._plot_model(lsq_fit_vals, model, ax1)

            # update values
            for key in lsq_fit_vals.keys():
                if key in self._full_params.keys():
                    self._full_params[key].value = lsq_fit_vals[key]

        # Evaluate and plot the model for the given set of parameters
        if apply_model and not lsq_fit:
            params = {key: val.value for key, val in self._full_params.items()}
            self._plot_model(params, model, ax1)

        fig.tight_layout()
        fig.show()

    def get_params_for_kriging(self) -> dict:
        """
        Collects all semi-variogram specific parameters
        needed for Kriging and returns them.

        Returns
        -------
        dict
            Dictionary with all parameters needed for Kriging.

        Notes
        -----
        This will return all parameters set by the ``get_widget``
        function, thus this function must be run first.
        """

        if self._model_drop:

            func, model_keys = self.get_model_info_from_str(self._model_drop.value)

            # select only those key, value pairs in the stored params that are needed by func
            model_params = dict((k, self._full_params[k].value) for k in model_keys if k in self._full_params.keys())

            return dict(s_v_model=func, s_v_params=model_params)

        else:
            raise RuntimeError("You must run get_widget() first!")

    def get_widget(self) -> ipywidgets.GridspecLayout:
        """
        Collects and returns a widget that allows one to
        plot the normalized semi-variogram values, apply
        a Least Squares fit to the normalized semi-variogram
        values (for a given model), and interactively plot a
        model for a given set of parameters.

        Returns
        -------
        ipywidgets.GridspecLayout
            An interactive layout for plotting and fitting.

        Notes
        -----
        If one would like a dictionary of all semi-variogram
        parameters needed for Kriging, one should consider
        using ``get_params_for_kriging``, after this function.
        """

        if not isinstance(self.gamma_normalized, np.ndarray):
            raise RuntimeError("You must calculate the semi-variogram first!")

        self._create_widgets()

        # construct interactive elements
        inter_out_exp = widgets.interactive_output(self._manage_plots,
                                                   {'model': self._model_drop,
                                                    'lsq_fit': self._lsq_toggle,
                                                    'apply_model': self._apply_mod_toggle
                                                    })
        inter_out_exp.observe(self._sill_box, 'value')

        # create grid layout for interactive display
        grid = widgets.GridspecLayout(8, 2, width='80%', height='auto')
        grid[0, 0] = self._model_drop
        grid[1, 0] = self._sill_box
        grid[2, 0] = self._len_scale_box
        grid[3, 0] = self._exp_pow_box
        grid[4, 0] = self._ls_hole_eff_box
        grid[5, 0] = self._nugget_box
        grid[6, 0] = self._apply_mod_toggle
        grid[7, 0] = self._lsq_toggle
        grid[:, 1:] = inter_out_exp

        return grid


