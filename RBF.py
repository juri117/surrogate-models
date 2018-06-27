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
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize


class RBF:

    def __init__(self, known_in, known_val):
        self._knownIn = np.array(known_in)
        self._knownVal = np.array(known_val)
        if len(self._knownIn.shape) == 1:
            self._knownIn = self._knownIn.reshape((1,self._knownIn.shape[0]))
        self._k = self._knownIn.shape[0]
        self._n = self._knownIn.shape[1]
        #else:
        #    self._k = 1
        #    self._n = self._knownIn.shape[0]

        self._coeff = None
        self._rbfConst = 1.
        self._rbf = gausRBF

    def update_param(self, rbf_const, rbf_name):
        self._rbfConst = rbf_const
        if rbf_name == 'gaus':
            self._rbf = gausRBF
        elif rbf_name == 'inverseBla':
            self._rbf = None
        else:
            print('WARNING: unknown rbf_name, I will just use gaus for you.')
            self._rbf = gausRBF
        self._calc_coefficiants()

    def _calc_coefficiants(self):
        paramCount, knownCount = np.array(self._knownIn).shape
        mat = np.zeros((knownCount, knownCount))
        radialMat = np.zeros((knownCount, knownCount))
        for iRow in range(0, knownCount):
            for iColumn in range(iRow, knownCount):
                radSum = 0.
                for iParam in range(0, paramCount):
                    radSum += (self._knownIn[iParam][iRow] - self._knownIn[iParam][iColumn]) ** 2.
                radius = math.sqrt(radSum)
                radialMat[iRow][iColumn] = radius
                mat[iRow][iColumn] = self._rbf(self._rbfConst, radius)
                mat[iColumn][iRow] = mat[iRow][iColumn]
        # print('radialMat:')
        # printMat(radialMat)
        # print('mat:')
        # printMat(mat)
        self._coeff = np.linalg.solve(mat, self._knownVal)
        return self._coeff

    def predict(self, x_pred):
        res = 0.
        for i in range(0, self._n):
            radSum = 0.
            for iParam in range(0, self._k):
                radSum += (x_pred[iParam] - self._knownIn[iParam][i]) ** 2.
            radius = math.sqrt(radSum)
            res += self._coeff[i] * self._rbf(self._rbfConst, radius)
        return res

def gausRBF(a, r):
    # return r
    # return math.e**(-((a*r)**2))
    return math.sqrt(1 + (a * r) ** 2)

# todo: other RBFunctions here...
