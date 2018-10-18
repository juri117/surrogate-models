__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional RadialBasisFunction
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
import math


class RBF:

    def __init__(self, known_in, known_val):
        self._known_in = np.array(known_in)
        self._known_val = np.array(known_val)
        if len(self._known_in.shape) == 1:
            self._known_in = self._known_in.reshape((self._known_in.shape[0], 1))
        self._k = self._known_in.shape[1]
        self._n = self._known_in.shape[0]
        self._coeff = None
        self._rbf_const = 1.
        self._rbf = gaus_rbf

    def train(self):
        pass
        #self.update_param(self._rbf_const, self._rbf)

    def update_param(self, rbf_const, rbf_name):
        self._rbf_const = rbf_const
        if rbf_name == 'lin':
            self._rbf = lin_rbf
        elif rbf_name == 'cubic':
            self._rbf = cubic_rbf
        elif rbf_name == 'gaus':
            self._rbf = gaus_rbf
        elif rbf_name == 'multi-quadratic':
            self._rbf = multi_quad_rbf
        elif rbf_name == 'inverse-multi-quadratic' or rbf_name == 'imq':
            self._rbf = inv_multi_quad_rbf
        else:
            print('WARNING: unknown rbf_name (' + rbf_name + '), I will just use gaus for you.')
            print('next time chose one of ["gaus", "multi-quadratic", "inverse-multi-quadratic"]')
            self._rbf = gaus_rbf
        self._calc_coefficiants()

    def _calc_coefficiants(self):
        mat = np.zeros((self._n, self._n))
        radial_mat = np.zeros((self._n, self._n))
        for i in range(0, self._n):
            for j in range(i, self._n):
                radSum = 0.
                for ik in range(0, self._k):
                    radSum += (self._known_in[i][ik] - self._known_in[j][ik]) ** 2.
                radius = math.sqrt(radSum)
                radial_mat[i][j] = radius
                mat[i][j] = self._rbf(self._rbf_const, radius)
                mat[j][i] = mat[i][j]
        self._coeff = np.linalg.solve(mat, self._known_val)
        return self._coeff

    def get_coeff(self):
        return self._coeff

    def predict(self, x_pred):
        res = 0.
        for i in range(0, self._n):
            rad_sum = 0.
            for ik in range(0, self._k):
                rad_sum += (x_pred[ik] - self._known_in[i][ik]) ** 2.
            radius = math.sqrt(rad_sum)
            res += self._coeff[i] * self._rbf(self._rbf_const, radius)
        return res

def lin_rbf(a, r):
    return r

def cubic_rbf(a, r):
    return r**3

def gaus_rbf(a, r):
    return math.e**(-((a*r)**2))

def multi_quad_rbf(a, r):
    return math.sqrt(1 + (a * r) ** 2)

def inv_multi_quad_rbf(a, r):
    return (1+r**2)**(a/2)
