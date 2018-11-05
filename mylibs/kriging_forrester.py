__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Kriging (alternativ implementation using non-numpy-methods)
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
import math
from scipy.optimize import minimize


class KrigingForrester:

    def __init__(self, known_in, known_val):
        """
        :param known_in: list of lists with input sample points
        :param known_val: list of results for the known_in
        """
        self._known_in = known_in
        self._known_val = known_val
        self._k = self._known_in.shape[0]
        self._n = self._known_in.shape[1]
        self._theta = 1. * np.ones((self._k, 1)).flatten()
        self._p = 2. * np.ones((self._k, 1)).flatten()
        self._cor_mat = None
        self._mat_u = None
        self._mu = None

    def update_param(self, theta, p):
        """
        updates the parameters of the surrogate model
        :param theta: vector of theta parameters for each entry one (from range [1e-5 .. 1e+5])
        :param p: vector of p parameters for each entry one (from range[1 ..2])
        :return: None
        """
        theta = np.array(theta)
        p = np.array(p)
        if theta.shape[0] != self._k or p.shape[0] != self._k:
            raise ValueError('wrong theta or p dimension, should be ' + str(self._k) + 'x1')
        self._calc_cormat()
        self._calc_mu()

    def _calc_cormat(self):
        cor_mat = np.zeros((self._n, self._n))
        for row in range(0, self._n):
            for column in range(row, self._n):
                sum = 0.
                for ik in range(0, self._k):
                    sum += self._theta[ik] * (abs(self._known_in[ik][row] - self._known_in[ik][column]) ** self._p[ik])
                cor_mat[row][column] = math.exp(-sum)
                cor_mat[column][row] = cor_mat[row][column]
        self._cor_mat = cor_mat
        self._mat_u = np.transpose(np.linalg.cholesky(self._cor_mat))
        return cor_mat

    def _calc_mu(self):
        one = np.ones((self._n, 1)).flatten()
        mu_top = np.transpose(one) @ (np.linalg.solve(self._mat_u, np.linalg.solve(np.transpose(self._mat_u), self._known_val)))
        mu_but = np.transpose(one) @ (np.linalg.solve(self._mat_u, np.linalg.solve(np.transpose(self._mat_u), one)))
        self._mu = mu_top / mu_but
        return self._mu

    def calc_likelihood(self):
        one = np.ones((self._n, 1)).flatten()
        ln_det_cor_mat = 2 * sum(np.log(abs(np.diag(self._mat_u))))
        sigma_sqr = (np.transpose(self._known_val - (one * self._mu)) @
                     np.linalg.solve(self._mat_u, np.linalg.solve(np.transpose(self._mat_u), (self._known_val - (one * self._mu))))) / self._n;
        neg_ln_like = (-1) * (-(self._n / 2) * np.log(sigma_sqr) - 0.5 * ln_det_cor_mat)
        if np.isnan(neg_ln_like):
            print('Error: nan')
        return neg_ln_like

    def _calc_likelihood_opti(self, params, *args):
        self.update_param(params, args[0])
        neg_ln_like = self.calc_likelihood()
        print(str(neg_ln_like))
        return neg_ln_like

    def optimize(self):
        x0 = [1., 1.]
        bnds = [(0.0001, 1000.), (0.0001, 1000.)]
        opt = {}
        opt['disp'] = True
        opt['maxiter'] = 99999
        res = minimize(self._calc_likelihood_opti, x0, args=(self._p), method='SLSQP', tol=1e-6, options=opt, bounds=bnds)
        self._theta = res.x

    def predict(self, x_pred):
        """
        predicts a value from the surrogate model
        :param x_pred: vector of input values
        :return: result value
        """
        one = np.ones((self._n, 1)).flatten()
        psi = np.ones((self._n, 1)).flatten()
        for i in range(0, self._n):
            sum = 0.
            for ik in range(0, self._k):
                sum += self._theta[ik] * (abs(self._known_in[ik][i] - x_pred[ik]) ** self._p[ik])
            psi[i] = math.exp(-sum)
        fx = self._mu + np.transpose(psi) @ \
             np.linalg.solve(self._mat_u,
                             np.linalg.solve(np.transpose(self._mat_u),
                                             (self._known_val - (one * self._mu))))
        return fx