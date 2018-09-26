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
        self._knownIn = np.array(known_in)
        self._knownVal = np.array(known_val)
        if len(self._knownIn.shape) == 1:
            self._knownIn = self._knownIn.reshape((self._knownIn.shape[0], 1))
        self._k = self._knownIn.shape[1]
        self._n = self._knownIn.shape[0]
        self._coeff = None
        self._rbfConst = 1.
        self._rbf = gausRBF

    def train(self):
        pass
        #self.update_param(self._rbfConst, self._rbf)

    def update_param(self, rbf_const, rbf_name):
        self._rbfConst = rbf_const
        if rbf_name == 'lin':
            self._rbf = linRBF
        elif rbf_name == 'cubic':
            self._rbf = cubicRBF
        elif rbf_name == 'gaus':
            self._rbf = gausRBF
        elif rbf_name == 'multi-quadratic':
            self._rbf = multiQuadRBF
        elif rbf_name == 'inverse-multi-quadratic' or rbf_name == 'imq':
            self._rbf = invMultiQuadRBF
        else:
            print('WARNING: unknown rbf_name (' + rbf_name + '), I will just use gaus for you.')
            print('next time chose one of ["gaus", "multi-quadratic", "inverse-multi-quadratic"]')
            self._rbf = gausRBF
        self._calc_coefficiants()

    def _calc_coefficiants(self):
        mat = np.zeros((self._n, self._n))
        radialMat = np.zeros((self._n, self._n))
        for i in range(0, self._n):
            for j in range(i, self._n):
                radSum = 0.
                for ik in range(0, self._k):
                    radSum += (self._knownIn[i][ik] - self._knownIn[j][ik]) ** 2.
                radius = math.sqrt(radSum)
                radialMat[i][j] = radius
                mat[i][j] = self._rbf(self._rbfConst, radius)
                mat[j][i] = mat[i][j]
        self._coeff = np.linalg.solve(mat, self._knownVal)
        return self._coeff

    def get_coeff(self):
        return self._coeff

    def predict(self, x_pred):
        res = 0.
        for i in range(0, self._n):
            radSum = 0.
            for ik in range(0, self._k):
                radSum += (x_pred[ik] - self._knownIn[i][ik]) ** 2.
            radius = math.sqrt(radSum)
            res += self._coeff[i] * self._rbf(self._rbfConst, radius)
        return res

def linRBF(a, r):
    return r

def cubicRBF(a, r):
    return r**3

def gausRBF(a, r):
    return math.e**(-((a*r)**2))

def multiQuadRBF(a, r):
    return math.sqrt(1 + (a * r) ** 2)

def invMultiQuadRBF(a, r):
    return (1+r**2)**(a/2)
