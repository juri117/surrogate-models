
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize


class KrigingForrester:

    def __init__(self, known_in, known_val):
        self._knownIn = known_in
        self._knownVal = known_val
        self._k = self._knownIn.shape[0]
        self._n = self._knownIn.shape[1]
        self._theta = 1. * np.ones((self._k, 1)).flatten()
        self._p = 2. * np.ones((self._k, 1)).flatten()
        self._corMat = None
        self._matU = None
        self._mu = None

    def update_param(self, theta, p):
        theta = np.array(theta)
        p = np.array(p)
        if theta.shape[0] != self._k or p.shape[0] != self._k:
            raise ValueError('wrong theta or p dimension, should be ' + str(self._k) + 'x1')
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
        self._matU = np.transpose(np.linalg.cholesky(self._corMat))
        return corMat

    def _calc_mu(self):
        one = np.ones((self._n, 1)).flatten()
        muTop = np.transpose(one) @ (np.linalg.solve(self._matU, np.linalg.solve(np.transpose(self._matU), self._knownVal)))
        muBut = np.transpose(one) @ (np.linalg.solve(self._matU, np.linalg.solve(np.transpose(self._matU), one)))
        self._mu = muTop / muBut
        return self._mu

    def calc_likelihood(self):
        one = np.ones((self._n, 1)).flatten()
        lnDetCorMat = 2 * sum(np.log(abs(np.diag(self._matU))))
        sigmaSqr = (np.transpose(self._knownVal - (one * self._mu)) @
                    np.linalg.solve(self._matU, np.linalg.solve(np.transpose(self._matU), (self._knownVal - (one * self._mu))))) / self._n;
        negLnLike = (-1) * (-(self._n / 2) * np.log(sigmaSqr) - 0.5 * lnDetCorMat)
        if np.isnan(negLnLike):
            print('Error: nan')
        return negLnLike

    def _calc_likelihood_opti(self, params, *args):
        self.update_param(params, args[0])
        NegLnLike = self.calc_likelihood()
        print(str(NegLnLike))
        return NegLnLike

    def optimize(self):
        x0 = [1., 1.]
        bnds = [(0.0001, 1000.), (0.0001, 1000.)]
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
        fx = self._mu + np.transpose(psi) @ \
             np.linalg.solve(self._matU,
                             np.linalg.solve(np.transpose(self._matU),
                                             (self._knownVal - (one * self._mu))))
        return fx