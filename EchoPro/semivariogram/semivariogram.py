import numpy as np
from scipy import special
import ipywidgets as widgets
import matplotlib.pyplot as plt
import warnings
import sys


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
        self._lsq_toggle = None

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
        bins

        Returns
        -------

        """

        self._center_bins = center_bins # TODO: put bins in ascending order

        x_diff = np.subtract.outer(self.x, self.x)
        y_diff = np.subtract.outer(self.y, self.y)

        field_head, field_tail = np.meshgrid(self.field, self.field, indexing='ij')
        field_diff_sqrd = np.power(field_head - field_tail, 2)

        # find the distance between points
        dis = np.sqrt(np.power(x_diff, 2) + np.power(y_diff, 2))

        # obtain the upper triangular portion of dis
        dis_2 = np.triu(dis, k=1)  # TODO: we can put this in a sparse form can remove dis_2 != 0, if implemented

        self.gamma_standardized = []
        for i in range(len(self._center_bins) - 1):
            # get indices of distances that in the lag
            bin_add = (self._center_bins[i + 1] - self._center_bins[i])/2.0
            ind_in_lag = np.argwhere((self._center_bins[i] - bin_add <= dis_2)
                                     & (dis_2 < self._center_bins[i + 1] + bin_add) & (dis_2 != 0))

            # indices in the lag
            x_ind = ind_in_lag[:, 0]
            y_ind = ind_in_lag[:, 1]

            # calculate the semi-variogram value
            gamma = 0.5 * np.mean(field_diff_sqrd[x_ind, y_ind])

            # standardize gamma by the standard deviation of the head
            # multiplied by the standard deviation of the tail
            std_head = np.std(field_head[x_ind, y_ind])
            std_tail = np.std(field_tail[x_ind, y_ind])
            self.gamma_standardized.append(gamma / (std_head * std_tail))

        self.gamma_standardized = np.array(self.gamma_standardized)

    def generalized_exp_bessel(self, lag_vec, sill, ls, exp_pow, ls_hole_eff, nugget):

        func = (sill - nugget)*(1 - np.exp(-(lag_vec/ls)**exp_pow)*special.j0(ls_hole_eff*lag_vec)) + nugget

        return func

    def get_widgets(self):

        self._model_drop = widgets.Dropdown(options=['Exponential', 'Generalized exponential-Bessel'],
                                            value='Exponential',
                                            description='Variogram Model',
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
                                                icon='line-chart')

    def plot_data(self,
                  sill=1.0,
                  len_scale=1.0,
                  model='Generalized exponential-Bessel',
                  lsq_fit=False):

        plt.close()
        fig = plt.figure(1, figsize=(8, 4), dpi=100)
        fig.canvas.manager.set_window_title('')
        plt.close(1)
        ax = fig.add_subplot()
        ax.plot(self._center_bins, self.gamma_standardized, 'bo')

        # if model == 'Generalized exponential-Bessel':
        #     ax.plot()
        # else:
        #     print("Not implemented!")

        if lsq_fit:
            ax.plot(0.002, 0.5, 'r.')

        fig.tight_layout()
        fig.show()

    def view_semi_variogram(self):

        if self._center_bins is None:
            print("You must calculate the semi-variogram first!")
            sys.exit()

        self.get_widgets()

        self._model_drop = None
        self._sill_box = None
        self._len_scale_box = None
        self._exp_pow_box = None
        self._ls_hole_eff_box = None
        self._nugget_box = None
        self._lsq_toggle = None

        inter_out_exp = widgets.interactive_output(self.plot_data,
                                                   {'sill': self._sill_box,
                                                    'len_scale': self._len_scale_box,
                                                    'model': self._model_drop,
                                                    'lsq_fit': self._lsq_toggle})

        grid = widgets.GridspecLayout(4, 2, width='80%', height='auto')

        grid[0, 0] = self._sill_box
        grid[1, 0] = self._len_scale_box
        grid[2, 0] = self._model_drop
        grid[3, 0] = self._lsq_toggle
        grid[:, 1:] = inter_out_exp

        return grid


