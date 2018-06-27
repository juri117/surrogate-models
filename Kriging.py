__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Kriging
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


class Kriging:

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

        self._theta = 1. * np.ones((self._k, 1)).flatten()
        self._p = 2. * np.ones((self._k, 1)).flatten()
        self._corMat = None
        self._coreMatInv = None
        self._mu = None

    def update_param(self, theta, p):
        self._theta = np.array(theta)
        self._p = np.array(p)
        #if theta.shape[0] != self._k or p.shape[0] != self._k:
        #    raise ValueError('wrong theta or p dimension, should be ' + str(self._k) + 'x1')
        self._calc_cormat()
        self._calc_mu()

    #calcs the correlation matrix
    def _calc_cormat(self):
        corMat = np.zeros((self._n, self._n))
        for row in range(0, self._n):
            for column in range(row, self._n):
                sum = 0.
                for ik in range(0, self._k):
                    sum += self._theta[ik] * (abs(self._knownIn[ik][row] - self._knownIn[ik][column]) ** self._p[ik])
                corMat[row][column] = math.exp(-sum)
                corMat[column][row] = corMat[row][column]
        self._corMat = corMat
        self._coreMatInv = np.linalg.inv(self._corMat)
        return corMat

    def _calc_mu(self):
        one = np.ones((self._n, 1)).flatten()
        self._mu = (np.transpose(one) @ self._coreMatInv @ self._knownVal) / (
                np.transpose(one) @ self._coreMatInv @ one)
        return self._mu

    #calcs the likelyhood
    def calc_likelihood(self):
        #lnDetCorMat = np.log(np.linalg.det(corMat))
        lnDetCorMat = np.linalg.slogdet(self._corMat)[1]
        one = np.ones((self._n, 1)).flatten()
        sigmaSqr = (np.transpose(self._knownVal - one * self._mu) @ self._coreMatInv @ (self._knownVal - one * self._mu)) / self._n
        negLnLike = (-1) * (-(self._n / 2) * np.log(sigmaSqr) - 0.5 * lnDetCorMat)
        if negLnLike == float('nan'):
            print('Error: nan')
        return negLnLike

    def _calc_likelihood_opti(self, params, *args):
        self.update_param(params, args[0])
        NegLnLike = self.calc_likelihood()
        print(str(NegLnLike))
        return NegLnLike

    def optimize(self):
        x0 = np.ones((self._k,1)).flatten()
        bnds = []
        for i in range(0, self._k):
            bnds.append((0.0001, 1000.))
        opt = {}
        opt['disp'] = True
        opt['maxiter'] = 99999
        res = minimize(self._calc_likelihood_opti, x0, args=(self._p), method='SLSQP', tol=1e-6, options=opt, bounds=bnds)
        self._theta = res.x

    def predict(self, x_pred):
        one = np.ones((self._n, 1)).flatten()
        psi = np.ones((self._n, 1)).flatten()
        for i in range(0, self._n):
            sum = 0.
            for ik in range(0, self._k):
                sum += self._theta[ik] * (abs(self._knownIn[ik][i] - x_pred[ik]) ** self._p[ik])
            psi[i] = math.exp(-sum)
        fx = self._mu + np.transpose(psi) @ self._coreMatInv @ (self._knownVal - one * self._mu)
        return fx
