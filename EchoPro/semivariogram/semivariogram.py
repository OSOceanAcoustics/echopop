import numpy as np
from scipy import special
import ipywidgets as widgets
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import inspect
import traitlets
import warnings
import sys

# default bounds for fitting the semi-variogram model
bnds_default = {'LB': {'sill': 0.0, 'ls': 0.0,
                       'exp_pow': 0.0, 'ls_hole_eff': 0.0, 'nugget': 0.0},
                'UB': {'sill': np.inf, 'ls': np.inf,
                       'exp_pow': np.inf, 'ls_hole_eff': 1e-13, 'nugget': 1e-13}}


class SemiVariogram:
    """
    This class contains a routine that calculates a
    standardized semi-variogram and routines for
    obtaining the best semi-variogram model for the
    estimated semi-variogram.

    Parameters
    ----------
    x : Numpy array
        A 1D array representing the x coordinate for
        the semi-variogram calculation.
    y : Numpy array
        A 1D array representing the y coordinate for
        the semi-variogram calculation.
    field : Numpy array
        A 1D array representing the field for the
        semi-variogram calculation (e.g. biomass density).
    """

    def __init__(self, x, y, field):

        self.x = x
        self.y = y
        self.field = field

        # bins provided to calculate_semi_variogram
        self._center_bins = None

        # standardized semi-variogram values
        self.gamma_standardized = None

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

    def calculate_semi_variogram(self, center_bins):
        """
        Calculates the semi-variogram standardized by the
        standard deviation of the head multiplied by
        the standard deviation of the tail for each
        lag. This calculation assumes that the mesh
        points (x,y) are isotropic and the search
        area is omnidirectional.

        Parameters
        ----------
        center_bins: Numpy Array  # TODO: fill in

        Returns
        -------

        """

        self._center_bins = np.sort(center_bins)

        # get upper triangular indices
        i_up_ind, j_up_ind = np.triu_indices(len(self.x), k=1)

        x_diff = np.subtract.outer(self.x, self.x)[i_up_ind, j_up_ind]
        y_diff = np.subtract.outer(self.y, self.y)[i_up_ind, j_up_ind]

        field_rep = np.tile(self.field, (len(self.field), 1))
        field_head = field_rep[j_up_ind, i_up_ind]
        field_tail = field_rep[i_up_ind, j_up_ind]
        field_diff = field_head - field_tail
        field_diff_sqrd = field_diff*field_diff

        # find the distance between points
        dis = np.sqrt(x_diff*x_diff + y_diff*y_diff)

        self.gamma_standardized = []
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

            # standardize gamma by the standard deviation of the head
            # multiplied by the standard deviation of the tail
            std_head = np.std(field_head[ind_in_lag])
            std_tail = np.std(field_tail[ind_in_lag])
            self.gamma_standardized.append(gamma / (std_head * std_tail))

        self.gamma_standardized = np.array(self.gamma_standardized)

    def get_widgets(self):

        self._model_drop = widgets.Dropdown(options=['Exponential', 'Generalized exponential-Bessel'],
                                            value='Generalized exponential-Bessel',
                                            description='Semi-variogram model',
                                            disabled=False, style={'description_width': 'initial'})

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

        self._lsq_toggle = widgets.ToggleButton(value=False,
                                                description='Run Least-squares fit',
                                                disabled=False, button_style='info',
                                                tooltip='Description',
                                                style={'description_width': 'initial'})

        self._apply_mod_toggle = widgets.ToggleButton(value=False,
                                                      description='Apply model',
                                                      disabled=False, button_style='info',
                                                      tooltip='Description',
                                                      icon='line-chart',
                                                      style={'description_width': 'initial'})

        self._full_params = {'sill': self._sill_box, 'ls': self._len_scale_box,
                             'exp_pow': self._exp_pow_box, 'ls_hole_eff': self._ls_hole_eff_box,
                             'nugget': self._nugget_box}

    @staticmethod
    def generalized_exp_bessel(lag_vec, sill, ls, exp_pow, ls_hole_eff, nugget):

        return (sill - nugget)*(1.0 - np.exp(-(lag_vec/ls)**exp_pow)*special.j0(ls_hole_eff*lag_vec)) + nugget

    @staticmethod
    def exponential(lag_vec, sill, ls, nugget):

        return (sill - nugget)*(1.0 - np.exp(-(lag_vec/ls))) + nugget

    def evaluate_model(self, params, model):

        if model == 'Generalized exponential-Bessel':
            func = self.generalized_exp_bessel
        elif model == 'Exponential':
            func = self.exponential

        model_keys = tuple(inspect.signature(func).parameters)
        sub_dict = dict((k, params[k]) for k in model_keys if k in params)
        return func(*tuple(sub_dict.values()))

    def apply_lsq_fit(self, model):

        if model == 'Generalized exponential-Bessel':
            func = self.generalized_exp_bessel
        elif model == 'Exponential':
            func = self.exponential

        model_keys = tuple(inspect.signature(func).parameters)

        low_bnds = [bnds_default['LB'][key] for key in model_keys if key in bnds_default['LB']]
        up_bnds = [bnds_default['UB'][key] for key in model_keys if key in bnds_default['UB']]

        fit_params, _ = curve_fit(func, self._center_bins, self.gamma_standardized,
                                  bounds=(low_bnds, up_bnds))

        return dict(zip(model_keys[1:], fit_params))

    def plot_model(self, all_inputs, model, ax):

        all_inputs['lag_vec'] = np.linspace(np.min(self._center_bins),
                                            np.max(self._center_bins), 100)
        model_vals = self.evaluate_model(all_inputs, model)

        ax.plot(all_inputs['lag_vec'], model_vals, 'r-')

        return ax

    def plot_data(self, model, lsq_fit, apply_model):

        plt.close(1)
        fig = plt.figure(1, figsize=(8, 4), dpi=100)
        fig.canvas.manager.set_window_title('')
        plt.close(1)
        ax1 = fig.add_subplot()

        ax1.plot(self._center_bins, self.gamma_standardized, 'bo')
        ax1.set_ylabel("gamma standardized")
        ax1.set_xlabel("lag")

        if lsq_fit:
            lsq_fit_vals = self.apply_lsq_fit(model)
            self.plot_model(lsq_fit_vals, model, ax1)

            # update values
            for key in lsq_fit_vals.keys():
                if key in self._full_params.keys():
                    self._full_params[key].value = lsq_fit_vals[key]

        if apply_model and not lsq_fit:
            params = {key: val.value for key, val in self._full_params.items()}
            self.plot_model(params, model, ax1)

        fig.tight_layout()
        fig.show()

    def view_semi_variogram(self):

        if self._center_bins is None:
            print("You must calculate the semi-variogram first!")
            sys.exit()

        self.get_widgets()

        inter_out_exp = widgets.interactive_output(self.plot_data,
                                                   {'model': self._model_drop,
                                                    'lsq_fit': self._lsq_toggle,
                                                    'apply_model': self._apply_mod_toggle
                                                    })

        inter_out_exp.observe(self._sill_box, 'value')

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


