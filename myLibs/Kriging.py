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

from utils.PlotHelper import PlotHelper
from utils.TimeTrack import TimeTrack

import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize
from scipy.optimize import basinhopping
from scipy import optimize


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

        try:
            self._coreMatInv = np.linalg.inv(corMat)
            self._corMat = corMat
        except:
            print('ERROR: could not calc np.linalg.inv')
        return corMat

    def _calc_mu(self):
        one = np.ones((self._n, 1)).flatten()
        self._mu = (np.transpose(one) @ self._coreMatInv @ self._knownVal) / (
                np.transpose(one) @ self._coreMatInv @ one)
        return self._mu

    #calcs the likelyhood
    def calc_likelihood(self):
        #lnDetCorMat = np.log(np.linalg.det(corMat))
        #ToDo: sometimes this throws an warning (C:\python\Python365_x64\lib\site-packages\numpy\linalg\linalg.py:1817: RuntimeWarning: invalid value encountered in slogdet sign, logdet = _umath_linalg.slogdet(a, signature=signature))
        lnDetCorMat = np.linalg.slogdet(self._corMat)[1]
        one = np.ones((self._n, 1)).flatten()
        #ToDo: sometimes this throws a warning (RuntimeWarning: invalid value encountered in log negLnLike = (-1) * (-(self._n / 2) * np.log(sigmaSqr) - 0.5 * lnDetCorMat))
        sigmaSqr = (np.transpose(self._knownVal - one * self._mu) @ self._coreMatInv @ (self._knownVal - one * self._mu)) / self._n
        if sigmaSqr < 0.:
            print('Error: neg sigmaSqr')
            return float('inf')
        negLnLike = (-1) * (-(self._n / 2) * np.log(sigmaSqr) - 0.5 * lnDetCorMat)
        if negLnLike == float('nan'):
            print('Error: nan')
        return negLnLike

    def _calc_likelihood_opti_theta_only(self, params, *args):
        self.update_param(params, args[0])
        NegLnLike = self.calc_likelihood()
        print(str(NegLnLike))
        return NegLnLike

    def _calc_likelihood_opti(self, params, *args):
        self.update_param(params[0:self._k], params[self._k:])
        NegLnLike = self.calc_likelihood()
        #print(str(NegLnLike))
        return NegLnLike

    def optimize_theta_only(self):
        x0 = np.ones((self._k,1)).flatten()
        bnds = []
        for i in range(0, self._k):
            bnds.append((0.0001, 1000.))
        opt = {}
        opt['disp'] = True
        opt['maxiter'] = 99999
        res = minimize(self._calc_likelihood_opti_theta_only, x0, args=self._p, method='SLSQP', tol=1e-6, options=opt, bounds=bnds)
        self._theta = res.x

    def optimize(self, init_guess=None):
        if init_guess is None:
            init_guess = [self._theta[0], self._theta[1], self._p[0], self._p[1]]
        bnds = []
        for i in range(0, self._k):
            bnds.append((1e-5, 1e+10))
        for i in range(0, self._k):
            bnds.append((1., 2.))

        tim1 = TimeTrack('LikelihoodOpti')
        tim1.tic()
        # SLSQP: proplem; find local min not glob. depending on init-vals
        #res = minimize(self._calc_likelihood_opti, init_guess, method='SLSQP', tol=1e-6, options={'disp': True, 'maxiter': 99999}, bounds=bnds)

        # random: not good enough... space is too big, cand find anything
        #res = optimize.differential_evolution(self._calc_likelihood_opti, bnds, maxiter=int(1e8))

        # basinhopping:
        bounds = BasinHoppingBounds()
        step = BasinHoppingStep()
        #L-BFGS-B
        minimizer_kwargs = dict(method='SLSQP', bounds=bnds, options={'disp': False, 'maxiter': 1e6}, tol=1e-4)
        res = basinhopping(self._calc_likelihood_opti,
                           init_guess,
                           minimizer_kwargs=minimizer_kwargs,
                           accept_test=bounds,
                           take_step=step,
                           niter=1000,
                           niter_success=100)
        self._theta = res.x[0:self._k]
        self._p = res.x[self._k:]
        tim1.toc()

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


    def plot_theta_likelihood_R2(self, ax=None):
        if self._k != 2:
            print('ERROR: plot_theta_likelihood_R2 only works with exactly 2 inputs')
            return

        opt_theta = self._theta
        thetas = np.logspace(-5, 9, num=50)
        likely_thet = np.zeros((len(thetas), len(thetas)))
        for i1 in range(0, len(thetas)):
            for i2 in range(0, len(thetas)):
                self.update_param([thetas[i1], thetas[i2]], self._p)
                likely_thet[i2][i1] = self.calc_likelihood()
        # restore original thetas
        self._theta = opt_theta
        self.update_param(self._theta, self._p)
        # plot it
        plt_theta = PlotHelper([r'$\theta_{1}$', r'$\theta_{2}$'], fancy=True, font_size=18, ax=ax)
        plt_theta.ax.set_xscale('log')
        plt_theta.ax.set_yscale('log')
        pcol = plt_theta.ax.pcolor(thetas, thetas, likely_thet, cmap='YlOrRd_r')
        #cbar = plt_theta.fig.colorbar(pcol)
        #cbar.set_label('neg. log. likelihood')
        plt_theta.ax.plot(self._theta[0], self._theta[1], 'wo', label='minimum')
        legend = plt_theta.finalize(width=6, height=5, legendLoc=4)
        legend.get_frame().set_facecolor('#000000')
        for text in legend.get_texts():
            text.set_color('#FFFFFF')
        #plt_theta.show()
        return pcol

    def plot_p_likelihood_R2(self, ax=None):
        if self._k != 2:
            print('ERROR: plot_p_likelihood_R2 only works with exactly 2 inputs')
            return
        opt_p = self._p
        ps = np.linspace(1., 2., num=50)
        likely_p = np.zeros((len(ps), len(ps)))
        for i1 in range(0, len(ps)):
            for i2 in range(0, len(ps)):
                self.update_param(self._theta, [ps[i1], ps[i2]])
                likely_p[i2][i1] = self.calc_likelihood()
        # restore original ps
        self._p = opt_p
        self.update_param(self._theta, self._p)
        # plot it
        plt_P = PlotHelper([r'$p_{1}$', r'$p_{2}$'], fancy=True, font_size=18, ax=ax)
        pcol = plt_P.ax.pcolor(ps, ps, likely_p, cmap='YlOrRd_r')
        #cbar = plt_P.fig.colorbar(pcol)
        #cbar.set_label('neg. log. likelihood')
        plt_P.ax.plot(self._p[0], self._p[1], 'wo', label='minimum')
        legend = plt_P.finalize(width=6, height=5, legendLoc=4)
        legend.get_frame().set_facecolor('#000000')
        for text in legend.get_texts():
            text.set_color('#FFFFFF')
        # plt_theta.show()
        return pcol

    def plot_likelihoods(self):
        figLike = plt.figure(figsize=(6, 4))
        ax1 = figLike.add_subplot(211)
        ax2 = figLike.add_subplot(212)

        pcol1 = self.plot_p_likelihood_R2(ax=ax1)
        pcol2 = self.plot_theta_likelihood_R2(ax=ax2)

        figLike.set_size_inches(6, 9)
        plt.tight_layout()

        # cbar = figLike.colorbar(pcol1)
        figLike.subplots_adjust(right=0.75)
        cbar_ax = figLike.add_axes([0.80, 0.15, 0.05, 0.78])
        figLike.colorbar(pcol2, cax=cbar_ax)
        # cbar_ax.set_label('neg. log. likelihood')
        figLike.text(0.93, 0.6, 'neg. log. likelihood', size=20, rotation=90.)


class BasinHoppingBounds(object):

    def __init__(self, xmax=[1e+10, 1e+10, 2., 2.], xmin=[1e-5, 1e-5, 1., 1.] ):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


class BasinHoppingStep(object):

    def __init__(self, stepsize=1.):
        self.stepsize = stepsize

    def __call__(self, x):
        s = self.stepsize
        x[0] = 10**np.random.uniform(-5, 10)
        x[1] = 10**np.random.uniform(-5, 10)
        x[2] = np.random.uniform(1., 2)
        x[3] = np.random.uniform(1., 2)
        #print('STEP: {:f}, {:f}, {:f}, {:f}'.format(x[0], x[1], x[2], x[3]))
        return x