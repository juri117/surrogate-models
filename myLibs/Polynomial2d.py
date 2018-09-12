__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Polynomial Model
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

from utils.PlotHelper import PlotHelper
from utils.TimeTrack import TimeTrack

import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy import optimize


class Polynomial2d:

    def __init__(self, known_in, known_val):
        self._knownIn = np.array(known_in)
        self._knownVal = np.array(known_val)
        if len(self._knownIn.shape) == 1:
            self._knownIn = self._knownIn.reshape((1, self._knownIn.shape[0]))
        self._k = self._knownIn.shape[0]
        self._n = self._knownIn.shape[1]
        self._order = 2
        if self._k != 1:
            raise('ERROR: Polynomial2d takes only 1 dimensional input')

    def update_param(self, order):
        self._order = order
        self._calc_vandermonde_mat()
        self._calc_weights()

    def _calc_vandermonde_mat(self):
        vander = np.zeros((self._n, self._order))
        for i in range(0, self._n):
            for o in range(0, self._order):
                vander[i][o] += self._knownIn[0][i]**o
        self._vander = vander

    def _calc_weights(self):
        # moore-penrose pseudo-inverse
        pin_vander = np.linalg.pinv(self._vander)
        weights = pin_vander @ self._knownVal
        self._weights = weights

    def predict(self, x_pred):
        fx = 0.
        iw = 0
        for o in range(0, self._order):
            fx += self._weights[iw] * x_pred[0]**o
            iw += 1
        return fx

    def get_order(self):
        return self._order

    def get_weights(self):
        return self._weights
