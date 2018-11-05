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

from myutils.plot_helper import PlotHelper
from myutils.time_track import TimeTrack
from mylibs.likeli_optimizer import LikeliOptimizer

import numpy as np
import math
from scipy.optimize import minimize
from scipy.optimize import basinhopping

VERBOSE = False


class Kriging:

    def __init__(self, known_in, known_val):
        """
        :param known_in: list of lists with input sample points
        :param known_val: list of results for the known_in
        """
        self._known_in = np.array(known_in)
        self._known_val = np.array(known_val)
        if len(self._known_in.shape) == 1:
            self._known_in = self._known_in.reshape((self._known_in.shape[0], 1))
        self._k = self._known_in.shape[1]
        self._n = self._known_in.shape[0]
        self._theta = 1. * np.ones((self._k, 1)).flatten()
        self._p = 2. * np.ones((self._k, 1)).flatten()
        self._cor_mat = None
        self._core_mat_inv = None
        self._mu = None
        self._mat_u = None

    def train(self):
        """
        trains the surrogate if available
        :return: None
        """
        self.optimize()

    def update_param(self, theta, p):
        """
        updates the parameters of the surrogate model
        :param theta: vector of theta parameters for each entry one (from range [1e-5 .. 1e+5])
        :param p: vector of p parameters for each entry one (from range[1 ..2])
        :return: None
        """
        self._theta = np.array(theta)
        self._p = np.array(p)
        self._calc_cormat()
        self._calc_mu()

    #calcs the correlation matrix
    def _calc_cormat(self):
        """
        :return: correlation matrix
        """
        cor_mat = np.zeros((self._n, self._n))
        for i in range(0, self._n):
            for j in range(i, self._n):
                sum = 0.
                for ik in range(0, self._k):
                    sum += self._theta[ik] * (abs(self._known_in[i][ik] - self._known_in[j][ik]) ** self._p[ik])
                cor_mat[i][j] = math.exp(-sum)
                cor_mat[j][i] = cor_mat[i][j]
        try:
            self._core_mat_inv = np.linalg.inv(cor_mat)
            self._cor_mat = cor_mat
        except Exception as e:
            if VERBOSE:
                print('ERROR: could not calc np.linalg.inv: ' + str(e))
        return cor_mat

    def _calc_mu(self):
        """
        :return: the factor mu
        """
        one = np.ones((self._n, 1)).flatten()
        self._mu = (np.transpose(one) @ self._core_mat_inv @ self._known_val) / (
                np.transpose(one) @ self._core_mat_inv @ one)
        return self._mu

    def calc_likelihood(self):
        """
        calculates the negative logarithmic likelihood
        :return: negative logarithmic likelihood (or infinity if an error appears)
        """
        ln_det_cor_mat = np.linalg.slogdet(self._cor_mat)[1]
        if np.isnan(ln_det_cor_mat):
            if VERBOSE:
                print('NaN Alarm')
            return float('inf')
        one = np.ones((self._n, 1)).flatten()
        sigma_sqr = (np.transpose(self._known_val - one * self._mu) @ self._core_mat_inv @ (self._known_val - one * self._mu)) / self._n
        if sigma_sqr < 0.:
            if VERBOSE:
                print('Error: neg sigmaSqr')
            return float('inf')
        neg_ln_like = (-1) * (-(self._n / 2) * np.log(sigma_sqr) - 0.5 * ln_det_cor_mat)
        if neg_ln_like == float('nan'):
            if VERBOSE:
                print('Error: nan')
            return float('inf')
        return neg_ln_like

    def _calc_likelihood_opti_theta_only(self, params, *args):
        self.update_param(params, args[0])
        neg_ln_like = self.calc_likelihood()
        return neg_ln_like

    def _calc_likelihood_opti(self, params, *args):
        self.update_param(params[0:self._k], params[self._k:])
        neg_ln_like = self.calc_likelihood()
        if self.records != None:
            self.records.append(params)
        return neg_ln_like

    def _calc_likelihood_opti_exp(self, params, *args):
        if np.isnan(params).any():
            return float('nan')
        exps = params[0:self._k]
        thetas = []
        for e in exps:
            thetas.append(10.**e)
        self.update_param(thetas, params[self._k:])
        neg_ln_like = self.calc_likelihood()
        if self.records != None:
            self.records.append(params)
        return neg_ln_like

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

    def optimize(self, init_guess=None, opti_algo='grid', record_data=False):
        """
        runs automatic optimization of thetas and ps
        :param init_guess: list of input values for an initial guess
        :param opti_algo: string for the algorithm to use 'grid' (self implemented LikeliOptimizer) or 'basin' (using scipy.optimize.basin-hopping)
        :param record_data: if True the test points of the optimizer gets recorded (this is needed for plot of optimizer path in plot_likelihoods)
        :return: None
        """
        timer = TimeTrack('optiTimer')
        self.records = None
        if record_data:
            self.records = []
        if 'basin' in opti_algo:
            # basinhopping:
            if init_guess is None:
                init_guess = []
                for t in self._theta:
                    init_guess.append(t)
                for p in self._p:
                    init_guess.append(p)
            bnds = []
            for i in range(0, self._k):
                bnds.append((-5., +5.))
            for i in range(0, self._k):
                bnds.append((1., 2.))
            bounds = BasinHoppingBounds(xmax=list(zip(*bnds))[1], xmin=list(zip(*bnds))[0])
            step = BasinHoppingStep()
            minimizer_kwargs = dict(method='SLSQP', bounds=bnds, options={'disp': False, 'maxiter': 5e3}, tol=1e-4)
            timer.tic()
            res = basinhopping(self._calc_likelihood_opti_exp,
                               init_guess,
                               minimizer_kwargs=minimizer_kwargs,
                               accept_test=bounds,
                               take_step=step,
                               niter=1000,
                               niter_success=100)
            timer.toc(print_it=True)
        elif 'grid' in opti_algo:
            skipper = LikeliOptimizer(debug=True)
            timer.tic()
            res = skipper.find(self._calc_likelihood_opti_exp, self._k)
            timer.toc(print_it=True)
        else:
            raise Exception('ERROR: unknown optimizer selected')
        exps = res.x[0:self._k]
        thetas = []
        for e in exps:
            thetas.append(10. ** e)
        if record_data:
            print('Kriging Likelihood optimization evaluations: {:d}'.format(len(self.records)))
        self.update_param(thetas, res.x[self._k:])

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
                sum += self._theta[ik] * (abs(self._known_in[i][ik] - x_pred[ik]) ** self._p[ik])
            psi[i] = math.exp(-sum)
        fx = self._mu + np.transpose(psi) @ self._core_mat_inv @ (self._known_val - one * self._mu)
        return fx

    def plot_theta_likelihood_r2(self, ax=None, pgf=False, opti_path=[]):
        """
        plot colormap of likelihood for theta1 and theta2
        :param ax: handle of the axis if this should be embedded in an existing plot
        :param pgf: store it as pgf file (for latex embedding)
        :param opti_path: show the path of the optimization as white crosses
        :return: the handle to the pcolor legend
        """
        if self._k != 2:
            print('ERROR: plot_theta_likelihood_R2 only works with exactly 2 inputs')
            return
        opt_theta = self._theta
        thetas = np.logspace(-5, 5, num=50)
        likely_thet = np.zeros((len(thetas), len(thetas)))
        for i1 in range(0, len(thetas)):
            for i2 in range(0, len(thetas)):
                self.update_param([thetas[i1], thetas[i2]], self._p)
                likely_thet[i2][i1] = self.calc_likelihood()
        # restore original thetas
        self._theta = opt_theta
        self.update_param(self._theta, self._p)
        # plot it
        plt_theta = PlotHelper([r'$\theta_{1}$', r'$\theta_{2}$'], fancy=False, ax=ax, pgf=pgf)
        plt_theta.ax.set_xscale('log')
        plt_theta.ax.set_yscale('log')
        pcol = plt_theta.ax.pcolor(thetas, thetas, likely_thet, cmap='YlOrRd_r')
        if len(opti_path) > 0:
            plt_theta.ax.plot(10**opti_path[:, 0], 10**opti_path[:, 1], '+', color='white', markeredgewidth=0.5, markersize=5, label='Optimierer-Pfad')
        plt_theta.ax.plot(self._theta[0], self._theta[1], 'x', color='black', label='Minimum', markersize=8, markeredgewidth=1.5)
        legend = plt_theta.finalize(width=6, height=5, legendLoc=4, show_legend=False)
        return pcol

    def plot_p_likelihood_r2(self, ax=None, pgf=False, opti_path=[]):
        """
        plot colormap of likelihood for p1 and p2
        :param ax: handle of the axis if this should be embedded in an existing plot
        :param pgf: store it as pgf file (for latex embedding)
        :param opti_path: show the path of the optimization as white crosses
        :return:
        """
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
        plt_P = PlotHelper(['$p_{1}$', '$p_{2}$'], fancy=False, ax=ax, pgf=pgf)
        pcol = plt_P.ax.pcolor(ps, ps, likely_p, cmap='YlOrRd_r')
        if len(opti_path) > 0:
            plt_P.ax.plot(opti_path[:, 0], opti_path[:, 1], '+', color='white', markeredgewidth=0.5, markersize=5, label='Optimierer-Pfad')
        plt_P.ax.plot(self._p[0], self._p[1], 'x', color='k', label='Minimum', markersize=8, markeredgewidth=1.5)
        legend = plt_P.finalize(width=6, height=5, legendLoc=4, show_legend=False)
        return pcol

    def plot_likelihoods(self, fancy=False, pgf=False, opti_path=[]):
        """
        creates a plot that shows the likelihood over p1, p2, theta1, theta2 (only supports two inputs)
        two colormaps represent the likelihoods for p and theta
        :param fancy: use latex in plot
        :param pgf: store it as pgf file (for latex embedding)
        :param opti_path: show the path of the optimization as white crosses
        :return: pointer to PlotHelper instance (if no errors)
        """
        if self._k != 2:
            print('ERROR: plot_p_likelihood_R2 only works with exactly 2 inputs')
            return
        path_p = []
        path_theta = []
        if len(opti_path) > 0:
            path_theta = np.array([opti_path[:, 0], opti_path[:, 1]]).T
            path_p = np.array([opti_path[:, 2], opti_path[:, 3]]).T

        plt_like = PlotHelper([], fancy=fancy, pgf=pgf)
        import matplotlib.pyplot as plt
        ax1 = plt_like.fig.add_subplot(121)
        ax2 = plt_like.fig.add_subplot(122)

        pcol1 = self.plot_p_likelihood_r2(ax=ax1, pgf=pgf, opti_path=path_p)
        pcol2 = self.plot_theta_likelihood_r2(ax=ax2, pgf=pgf, opti_path=path_theta)
        like_min = min(min(pcol1._A), min(pcol2._A))
        like_max = max(max(pcol1._A), max(pcol2._A))
        pcol1.set_clim(like_min, like_max)
        pcol2.set_clim(like_min, like_max)

        plt_like.fig.set_size_inches(6, 3)
        plt.tight_layout()

        plt_like.fig.subplots_adjust(right=0.85)
        plt_like.fig.subplots_adjust(bottom=0.3)
        cbar_ax = plt_like.fig.add_axes([0.88, 0.15, 0.02, 0.78])
        plt_like.fig.colorbar(pcol2, cax=cbar_ax)
        plt_like.fig.text(0.97, 0.7, 'neg. log. Likelihood', size=plt_like.FONT_SIZE, rotation=90.)
        handles, labels = ax1.get_legend_handles_labels()
        legend = plt_like.fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0.01), ncol=2, fancybox=True)
        legend.get_frame().set_facecolor('#A3A3A3')
        for text in legend.get_texts():
           text.set_color('#000000')
        return plt_like

    def get_p(self):
        return self._p

    def get_theta(self):
        return self._theta


"""
helper classes for basin hopping
"""
class BasinHoppingBounds(object):

    def __init__(self, xmax=[5., 5., 2., 2.], xmin=[-5., -5., 1., 1.]):
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
        for i in range(0, len(x)):
            if i < len(x) / 2:
                # theta
                x[i] = 10**np.random.uniform(-5, 5)
            else:
                # p
                x[i] = np.random.uniform(1., 2)
        return x