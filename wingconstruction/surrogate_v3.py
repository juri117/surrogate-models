__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :creates the surrogate for wing-structure test
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import sys
import os
from scipy.optimize import minimize
import numpy as np
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyKriging')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/inspyred')

from wingconstruction.wingutils.constants import Constants
from wingconstruction.multi_run import MultiRun
from myutils.time_track import TimeTrack
from myutils.plot_helper import PlotHelper
from mylibs.kriging import Kriging
from mylibs.rbf import RBF
from mylibs.polynomial import Polynomial
from mylibs.interface.rbf_scipy import RBFscipy
from mylibs.latin_hyper_cube import LatinHyperCube
from mylibs.halton import Halton
from mylibs.interface.opti_latin_hyper import OptiLatinHyper
from mylibs.structured_sample import StructuredSample
from mylibs.validation import Validation
from mylibs.validation import ValidationResults
from wingconstruction.fem.wing_construction_v4 import WingConstruction
from wingconstruction.wingutils.defines import *

FANCY_PLOT = True
RESULTS_FILE = '/2drun_2018-08-23_16_49_18_final01_cruiseLoad.csv'

'''
performs full surrogate analysis and comparison to real FEM-Model
:param sampling_type SAMPLE_LATIN, SAMPLE_HALTON, SAMPLE_STRUCTURE
:param sample_point_count amount of sample points
:param surro_type SURRO_KRIGING, SURRO_RBF
:param use_abaqus if False Calculix will be used
:param pgf if set to True LaTex ready plot files will be written to disk
:param show_plots weather plots will be displayed at runtime
:return SurroResults with all the results
'''
class Surrogate:

    def __init__(self, use_abaqus=False, pgf=False, show_plots=True, force_recalc=False, scale_it=True):
        self.use_abaqus = use_abaqus
        self.force_recalc = force_recalc
        self.pgf = pgf
        self.show_plots = show_plots
        self.scale_it = scale_it

        self.sampling_type = None
        self.surro_type = None
        self.ribs = None
        self.shell = None
        self.stress = None
        #self.weights = None
        self.ribs_s = None
        self.shell_s = None
        self.known_params = None
        self.known_params_s = None
        self.known_stress = None
        self.vali_params = None
        self.vali_params_s = None
        self.vali_values = None
        self.offset_rib = 0.
        self.offset_shell = 0.
        self.scale_rib = 1.
        self.scale_shell = 1.
        self.update_params = None
        self.surro_class = None

        self.multi = MultiRun(use_calcu=True, use_aba=True, non_liner=False, force_recalc=force_recalc)
        self.results = SurroResults()
        self.prepare(force_recalc)

    def auto_run(self, sampling_type, sample_point_count, surro_type, run_validation=True, auto_fit=True, params=[], sequential_runs=0):
        suc = self.generate_sampling_plan(sampling_type, sample_point_count)
        if not suc:
            return self.results, None
        self.run_validation_points()
        for i in range(0, sequential_runs + 1):
            if self.results.optimum_weight > 0.:
                self.known_params = np.append(self.known_params,
                                              [[self.results.optimum_rib, self.results.optimum_shell]], axis=0)
                self._generate_scaled_sampling_points()

            self.run_fem_calculation()

            #fit_time = TimeTrack('FitTime')
            #fit_time.tic()
            if surro_type == SURRO_POLYNOM and auto_fit:
                suc = self.auto_fit_poly()
            elif surro_type == SURRO_RBF and auto_fit:
                suc = self.auto_fit_rbf(params=params)
            else:
                suc = self.train_model(surro_type, params=params)
            #self.results.runtime = fit_time.toc()
            if not suc:
                return self.results, None
            self.optimize()

        if run_validation:
            self.run_validation_points(flip_points=True)
            self.run_validation(full_validation=True)

        self.print_results()
        if self.show_plots:
            self.plot_it()
        return self.results, self.surro

    def auto_fit_poly(self):
        orders = range(1, 15)
        rmse = np.zeros((len(orders)))
        for i in range(0, len(orders)):
            self.train_model(SURRO_POLYNOM, [orders[i]])
            self.run_validation(full_validation=False)
            rmse[i] = self.results.vali_results.rmse
        best_order = orders[np.argmin(rmse)]
        self.results.opti_params = [best_order]
        self.train_model(SURRO_POLYNOM, [best_order])
        return True

    def auto_fit_rbf(self, params=[]):
        rbf_func = 'gaus' # 'gaus' 'multi-quadratic'
        init_a = 1.
        if len(params) == 2:
            init_a = params[0]
            rbf_func = params[1]
        res = minimize(self._opti_rbf, init_a, args=rbf_func, method='SLSQP', tol=1e-8,
                       options={'disp': False, 'maxiter': 999}, bounds=[(0.01, 5.)])
        self.results.opti_params = [res.x[0], rbf_func]
        self.train_model(SURRO_RBF, [res.x[0], rbf_func])
        return True

    def _opti_rbf(self, a, rbf_func):
        self.train_model(SURRO_RBF, [a, rbf_func])
        self.run_validation(full_validation=False)
        return self.results.vali_results.rmse / max_shear_strength

    def prepare(self, force_recalc=False):
        ##################################################
        # collect data
        ribs, shell, stress, disp, weights = self.multi.read_data_file(RESULTS_FILE, use_abaqus=self.use_abaqus)
        # input 0 is ribs
        # input 1 is shell_thick
        # reduce data range
        i_rib_0 = 0
        i_rib_n = len(ribs) - 1
        for i in range(0, len(ribs)):
            if range_rib[0] >= ribs[i]:
                i_rib_0 = i
            if range_rib[1] <= ribs[i]:
                i_rib_n = i
                break
        i_shell_0 = 0
        i_shell_n = len(shell) - 1
        for i in range(0, len(shell)):
            if range_shell[0] >= shell[i]:
                i_shell_0 = i
            if range_shell[1] <= shell[i]:
                i_shell_n = i
                break
        new_rib = []
        new_shell = []
        new_stress = np.zeros((i_shell_n - i_shell_0 + 1, i_rib_n - i_rib_0 + 1))
        newWeight = np.zeros((i_shell_n - i_shell_0 + 1, i_rib_n - i_rib_0 + 1))
        for r in range(0, i_rib_n - i_rib_0 + 1):
            new_rib.append(ribs[i_rib_0 + r])
            for s in range(0, i_shell_n - i_shell_0 + 1):
                new_stress[s][r] = stress[i_shell_0 + s][i_rib_0 + r]
                newWeight[s][r] = weights[i_shell_0 + s][i_rib_0 + r]
        for s in range(0, i_shell_n - i_shell_0 + 1):
            new_shell.append(shell[i_shell_0 + s])
        self.ribs = np.array(new_rib)
        self.shell = np.array(new_shell)
        self.stress = new_stress
        weights = newWeight
        # calc offset and factors for scaling the inputs (some surrogates like scaled inputs -> RBF)
        if self.scale_it:
            self.offset_rib = range_rib[0]
            self.offset_shell = range_shell[0]
            self.scale_rib = range_rib[1] - range_rib[0]
            self.scale_shell = range_shell[1] - range_shell[0]
        self.ribs_s = (self.ribs - self.offset_rib) / self.scale_rib
        self.shell_s = (self.shell - self.offset_shell) / self.scale_shell

    def generate_sampling_plan(self, sampling_type, sample_point_count):
        self.sampling_type = sampling_type
        ##################################################
        # sample plan
        if sampling_type == SAMPLE_LATIN:
            sam = LatinHyperCube()
            sample_points = sam.generate_sample_plan(sample_point_count, 2, [range_rib, range_shell])
        elif sampling_type == SAMPLE_HALTON:
            sam = Halton()
            sample_points = sam.generate_sample_plan(sample_point_count, 2, [range_rib, range_shell], base=[2, 19])
        elif sampling_type == SAMPLE_STRUCTURE:
            sam = StructuredSample()
            sample_points = sam.generate_sample_plan(sample_point_count, 2, [range_rib, range_shell])
        elif sampling_type == SAMPLE_OPTI_LATIN_HYPER:
            sam = OptiLatinHyper()
            sample_points = sam.generate_sample_plan(sample_point_count, 2, [range_rib, range_shell])
        else:
            print('unknown sample plan selected')
            self.results.errorStr = 'unknown sample plan selected'
            return False
        # make the ribs be int
        for i in range(0, len(sample_points)):
            sample_points[i][0] = int(round(sample_points[i][0]))
        self.known_params = np.array(sample_points)
        self._generate_scaled_sampling_points()
        print('sample plan using {:d} known values'.format(len(self.known_params[:, 0])))
        return True

    def _generate_scaled_sampling_points(self):
        self.known_params_s = np.zeros(self.known_params.shape)
        self.known_params_s[:, 0] = (self.known_params[:, 0] - self.offset_rib) / self.scale_rib
        self.known_params_s[:, 1] = (self.known_params[:, 1] - self.offset_shell) / self.scale_shell

    def run_fem_calculation(self):
        ##################################################
        # FEM calculation, collecting results
        self.known_stress = self.multi.run_sample_points(self.known_params[:, 0], self.known_params[:, 1], use_abaqus=self.use_abaqus)

    def train_model(self, surro_type, params=[]):
        self.surro_type = surro_type
        ##################################################
        # build surrogate model and fit it
        # will be needed later for Validation
        # scale the inputs if needed
        if surro_type == SURRO_KRIGING:
            self.surro_class = Kriging
            self.surro = Kriging(self.known_params_s, self.known_stress)
            # prev stored results:
            self.surro.update_param([0.002261264770141511, 277826.21903867245], [1.8766170168043503, 1.9959876593551822])
            print('starting Likelihood optimization')
            opti_algo = 'grid'  # 'grid', 'basin'
            self.surro.optimize(opti_algo=opti_algo, record_data=True)
            if self.show_plots:
                pltLike = self.surro.plot_likelihoods(fancy=FANCY_PLOT, pgf=self.pgf, opti_path=np.array(self.surro.records))
                pltLike.save(Constants().PLOT_PATH + 'wingSurroLikely' + opti_algo + '.pdf')
                minLike = self.surro.calc_likelihood()
                #print('minLike = ' + str(minLike))
                #print('@theta1 = ' + str(self.surro.get_theta()[0]))
                #print('@theta2 = ' + str(self.surro.get_theta()[1]))
                #print('@p1 = ' + str(self.surro.get_p()[0]))
                #print('@p2 = ' + str(self.surro.get_p()[1]))
        elif surro_type == SURRO_RBF:
            self.surro_class = RBF
            self.surro = RBF(self.known_params_s, self.known_stress)
            a = 1.5
            rbf_func = 'gaus' # 'gaus' 'multi-quadratic'
            if params != []:
                a = params[0]
                if len(params) > 1:
                    rbf_func = params[1]
            self.surro.update_param(a, rbf_func)
            self.update_params = [a, rbf_func]
            #if self.show_plots:
            #    print('coeff1 = ' + str(self.surro.get_coeff()[0]))
            #    print('coeff2 = ' + str(self.surro.get_coeff()[1]))
        elif surro_type == SURRO_POLYNOM:
            self.surro_class = Polynomial
            self.surro = Polynomial(self.known_params_s, self.known_stress)
            o = 3
            if params != []:
                o = params[0]
            self.surro.update_param(o)
            #self.surro.generate_formula()
            self.update_params = [o]
        elif surro_type == SURRO_PYKRIGING:
            from pyKriging.krige import kriging as PyKriging
            self.surro_class = PyKriging
            self.surro = PyKriging(self.known_params_s, self.known_stress)
            self.surro.train()
        elif surro_type == SURRO_RBF_SCIPY:
            self.surro_class = RBFscipy
            self.surro = RBFscipy(self.known_params_s, self.known_stress)
            a = .8
            self.surro.update_param(a, 'linear')  # 'multi-quadratic')
            self.update_params = [a, 'linear']
        else:
            print('unknown surrogate type selected')
            self.results.errorStr = 'unknown surrogate type selected'
            return False
        return True

    def run_validation_points(self, flip_points=False):
        factor = 1
        if flip_points:
            factor = -1
        center_r = (range_rib[1] + range_rib[0]) * 0.5
        center_s = (range_shell[1] + range_shell[0]) * 0.5
        d_r = (range_rib[1] - range_rib[0]) / 2.
        d_s = (range_shell[1] - range_shell[0]) / 2.
        self.vali_params = np.array([[center_r, center_s],
                                [center_r - int(round(d_r / 4.) * factor), center_s - (factor * d_s / 4.)],
                                [center_r - int(round(d_r / 3.) * factor), center_s + (factor * d_s / 3.)],
                                [center_r + int(round(d_r / 2.7) * factor), center_s + (factor * d_s / 2.7)],
                                [center_r + int(round(d_r / 2.2) * factor), center_s - (factor * d_s / 2.2)]])
        self.vali_params_s = np.zeros(self.vali_params.shape)
        self.vali_params_s[:, 0] = (self.vali_params[:, 0] - self.offset_rib) / self.scale_rib
        self.vali_params_s[:, 1] = (self.vali_params[:, 1] - self.offset_shell) / self.scale_shell
        self.vali_values = self.multi.run_sample_points(self.vali_params.T[0], self.vali_params.T[1],
                                                   use_abaqus=self.use_abaqus)

    def run_validation(self, full_validation=False):
        ##################################################
        # validate
        vali = Validation()
        if full_validation:
            p_x, p_y = np.meshgrid(self.ribs, self.shell)
            params = np.array([p_x.flatten(), p_y.flatten()]).T
            params_s = np.zeros(params.shape)
            params_s[:, 0] = (params[:, 0] - self.offset_rib) / self.scale_rib
            params_s[:, 1] = (params[:, 1] - self.offset_shell) / self.scale_shell
            values = self.stress.flatten()
            vali_r = vali.run_full_analysis(params_s, values,
                                            self.known_params_s, self.known_stress,
                                            self.vali_params_s, self.vali_values,
                                            self.surro.predict, self.surro_class, update_params=self.update_params)
            self.results.vali_results = vali_r
        else:
            rmse = vali.calc_rmse(self.vali_params_s, self.vali_values, self.surro.predict)
            self.results.vali_results.rmse = rmse
        if self.show_plots and full_validation:
            deri_plot = PlotHelper(['Rippen', 'Blechdicke in mm'], fancy=FANCY_PLOT, pgf=self.pgf)
            dev = np.zeros(self.stress.shape)
            for xi in range(0, len(self.ribs_s)):
                for yi in range(0, len(self.shell_s)):
                    devi = (abs(self.stress[yi][xi] - self.surro.predict([self.ribs_s[xi], self.shell_s[yi]])) / np.array(
                        self.stress).mean()) * 100.
                    dev[yi][xi] = devi
            pcol = deri_plot.ax.pcolor(self.ribs, np.array(self.shell) * 1000, dev, cmap='YlOrRd', alpha=0.7)
            pcol.set_clim(0, 5.)
            cbar = deri_plot.fig.colorbar(pcol)
            deri_plot.ax.plot(self.known_params[:, 0], self.known_params[:, 1] * 1000, 'bo', label='Stützstellen')
            deri_plot.ax.plot(self.vali_params[:, 0], self.vali_params[:, 1] * 1000, 'o', color='fuchsia',
                              label='Vali.-Punkte')
            # deri_plot.ax.plot([opti_ribs[best_i]], [opti_shell[best_i] * 1000.], 'rx',
            #                  markersize=12, markeredgewidth=5, label='glob. Optimum')
            deri_plot.ax.invert_yaxis()
            deri_plot.finalize(width=6., height=4., legendLoc=8, legendNcol=3, bbox_to_anchor=(0.5, -0.38),
                               tighten_layout=True)
            deri_plot.save(
                Constants().PLOT_PATH + 'wingSurro_deri_{:s}_{:s}.pdf'.format(SAMPLE_NAMES[self.sampling_type],
                                                                                  SURRO_NAMES[self.surro_type]))

    @staticmethod
    def shell_predict(shell_thick, surro_inst, rib_num):
        stress_val = surro_inst.predict([rib_num, shell_thick])
        return stress_val - max_shear_strength

    def optimize(self):
        ##################################################
        # optimize
        opti_ribs = []
        opti_shell = []
        opti_stress = []
        opti_weights = []
        used_ribs = np.array(range(range_rib[0], range_rib[1] + 1))
        used_ribs_s = (used_ribs - self.offset_rib) / self.scale_rib
        for i in range(0, len(used_ribs_s)):
            # SLSQP: proplem; find local min not glob. depending on init-vals
            init_guess = (min(self.known_params_s[:, 1]) + max(self.known_params_s[:, 1])) / 2
            # bnds = [(min(known_shell), max(known_shell))]
            # res = minimize(shell_predict, init_guess, args=[krig, ribs[i]], method='SLSQP', tol=1e-6, options={'disp': True, 'maxiter': 99999}, bounds=bnds)
            # opti_shell.append(res.x[0])
            try:
                root_s = optimize.newton(self.shell_predict, init_guess, args=[self.surro, used_ribs_s[i]])
                root = (root_s * self.scale_shell) + self.offset_shell
                root_stress = self.surro.predict([used_ribs_s[i], root_s])
                if root_stress < max_shear_strength * 1.05: #this check is needed if the surrogate does not cross the max stress at all at this ribnumber
                    opti_ribs.append(used_ribs[i])
                    opti_shell.append(root)
                    opti_stress.append(root_stress)
                    weight = WingConstruction.calc_weight_stat(wing_length,
                                                               chord_length * 0.4,
                                                               chord_height,
                                                               used_ribs[i],
                                                               root,
                                                               density)
                    opti_weights.append(weight)
            except Exception as e:
                print(e)
        # exclude model edges from opti vals
        opti_ribs = opti_ribs[:-1]
        opti_ribs = opti_ribs[1:]
        opti_shell = opti_shell[:-1]
        opti_shell = opti_shell[1:]
        opti_stress = opti_stress[:-1]
        opti_stress = opti_stress[1:]
        opti_weights = opti_weights[:-1]
        opti_weights = opti_weights[1:]
        if len(opti_weights) > 0:
            best_i = opti_weights.index(min(opti_weights))
            if self.show_plots:
                optWeightPlot = PlotHelper(['ribs', 'weight'], pgf=self.pgf)
                optWeightPlot.ax.plot(opti_ribs, opti_weights, '-', color='dodgerblue')
                optWeightPlot.ax.plot([opti_ribs[best_i]], opti_weights[best_i], 'rx', label='minimum')
                import matplotlib.ticker as ticker
                optWeightPlot.ax.xaxis.set_major_locator(ticker.IndexLocator(base=2, offset=0))
                optWeightPlot.finalize(height=2)
            self.results.optimum_rib = opti_ribs[best_i]
            self.results.optimum_shell = opti_shell[best_i]
            self.results.optimum_weight = opti_weights[best_i]
            self.results.optimum_stress = opti_stress[best_i]
            self.results.opti_curve = [opti_ribs, opti_shell, opti_stress, opti_weights]

    def plot_it(self, display_plots=True):
        ##################################################
        # plot it
        # convert all to np arrays
        # shell = np.array(self.shell)
        # opti_shell = np.array(opti_shell)
        # known_shell = np.array(known_shell)

        if False:
            plot3dw = PlotHelper(['Rippen', 'Blechdicke in mm', 'Gewicht in kg'], fancy=False, pgf=pgf)
            plotX, plotY = np.meshgrid(self.ribs, self.shell * 1000)
            surf = plot3dw.ax.plot_wireframe(plotX,
                                             plotY,
                                             self.weights,
                                             color='blue',
                                             label='weight',
                                             rcount=20,
                                             ccount=20,
                                             linewidths=1,
                                             alpha=0.5)
            samplePoints = plot3dw.ax.plot(self.known_params[:, 0], self.known_params[:, 1] * 1000., 'bo', label='Stützstellen')
            plot3dw.finalize(show_legend=True)

        plot3d = PlotHelper([r'Rippen', r'Blechdicke in mm', r'Mises in Pa'], fancy=FANCY_PLOT, pgf=self.pgf)
        # plot FEM data as lines
        for i in range(0, len(self.ribs)):
            if i == 0:
                plot3d.ax.plot(np.ones((len(self.shell))) * self.ribs[i], self.shell * 1000., self.stress[:, i], 'g-', lw=3.,
                               label=u'FEM Daten')
            else:
                plot3d.ax.plot(np.ones((len(self.shell))) * self.ribs[i], self.shell * 1000., self.stress[:, i], 'g-', lw=3.)
        # plot surrogate model as wireframe
        ribs_sample = np.linspace(min(self.ribs_s), max(self.ribs_s), 200)
        shell_sample = np.linspace(min(self.shell_s), max(self.shell_s), 200)
        surro_short_name = SURRO_NAMES[self.surro_type][:3]
        if len(SURRO_NAMES[self.surro_type]) > 3:
            surro_short_name += '.'
        surro_plot = plot3d.plot_function_3d(self.surro.predict, ribs_sample, shell_sample,
                                             r'$\widehat{f}_{' + surro_short_name + '}$', color='b',
                                             scale=[self.scale_rib, self.scale_shell * 1000., 1.],
                                             offset=[self.offset_rib, self.offset_shell * 1000, 0.])
        samplePoints = plot3d.ax.plot(self.known_params[:, 0], self.known_params[:, 1] * 1000., self.known_stress, 'bo',
                                      label=u'Stützstellen')

        if self.results.opti_curve != []:
            # plot limit load line
            plot3d.ax.plot(self.results.opti_curve[0], np.array(self.results.opti_curve[1]) * 1000., self.results.opti_curve[2], 'k--', lw=3., label=u'max. Mises Stress')
            # plot optimal point
            plot3d.ax.plot([self.results.optimum_rib], [self.results.optimum_shell * 1000.], [self.results.optimum_stress], 'rx',
                       markersize=12,
                       markeredgewidth=5, label='glob. Optimum')
            plot3d.ax.plot([17], [2.565], [max_shear_strength], 'kx')
        plot3d.ax.locator_params(nbins=7, axis='x')
        plot3d.ax.locator_params(nbins=5, axis='y')

        plot3d.ax.set_zlim3d(np.min(np.array(self.stress)), max_shear_strength * 1.2)
        plot3d.ax.set_ylim3d(np.min(np.array(self.shell)) * 1000., np.max(np.array(self.shell)) * 1000.)

        plot3d.finalize(height=4, width=6, legendLoc=8, legendNcol=3, bbox_to_anchor=(0.5, -0.33))
        plot3d.ax.view_init(18, 105)
        plot3d.ax.invert_xaxis()
        # plot3d.ax.zaxis.offsetText.set_visible(True)
        # offset_z = plot3d.ax.zaxis.get_major_formatter().get_offset()
        offset_z = '1e8'
        plot3d.ax.set_zlabel('Mises in Pa x ' + offset_z, fontdict=plot3d.font, labelpad=plot3d.labelpad)
        plot3d.ax.zaxis.offsetText.set_visible(False)
        plot3d.save(Constants().PLOT_PATH + 'wingSurro_{:s}_{:s}.pdf'.format(SAMPLE_NAMES[self.sampling_type],
                                                                             SURRO_NAMES[self.surro_type]))
        if display_plots:
            plot3d.show()

    def print_results(self):
        print('optimum:')
        print('ribs: {:f}'.format(self.results.optimum_rib))
        print('shell: {:f}'.format(self.results.optimum_shell))
        print('weight: {:f}'.format(self.results.optimum_weight))
        print('stress: {:f}'.format(self.results.optimum_stress))
        print('opti params: ' + str(self.results.opti_params))
        print('runntime: {:f}'.format(self.results.runtime))
        print('---')
        print('deviation: {:f}'.format(100. * self.results.vali_results.deviation))
        print('rmse: {:f} ({:f}%)'.format(self.results.vali_results.rmse, 100. * self.results.vali_results.rmse / max_shear_strength))
        print('mae: {:f} ({:f}%)'.format(self.results.vali_results.mae, 100. * self.results.vali_results.mae / max_shear_strength))
        print('press: {:f} ({:f}%)'.format(self.results.vali_results.press, 100. * self.results.vali_results.press / max_shear_strength))

class SurroResults:

    def __init__(self):
        self.optimum_rib = 0.
        self.optimum_shell = 0.
        self.optimum_weight = 0.
        self.optimum_stress = 0.
        self.runtime = 0.
        self.errorStr = '-'
        self.vali_results = ValidationResults()
        self.opti_curve = []
        self.opti_params = []



#def surrogate_analysis(sampling_type, sample_point_count, surro_type, use_abaqus=False, pgf=False, show_plots=True, run_validation=True, scale_it=True):
#    return results, surro





if __name__ == '__main__':
    PGF = False
    SHOW_PLOT = True
    # SAMPLE_LATIN, SAMPLE_HALTON, SAMPLE_STRUCTURE, SAMPLE_OPTI_LATIN_HYPER
    # SURRO_KRIGING, SURRO_RBF, SURRO_POLYNOM, SURRO_PYKRIGING, SURRO_RBF_SCIPY
    if False:
        sur = Surrogate(use_abaqus=True, pgf=PGF, show_plots=SHOW_PLOT, scale_it=True)
        res, _ = sur.auto_run(SAMPLE_HALTON, 16, SURRO_KRIGING, run_validation=False, params=[1., 'gaus'], auto_fit=True, sequential_runs=0) # 'gaus' 'multi-quadratic'
    else:
        SHOW_PLOT = False
        SAMPLING = SAMPLE_HALTON
        VALIDATION = True
        POINTS = 17
        surP = Surrogate(use_abaqus=True, pgf=PGF, show_plots=SHOW_PLOT, scale_it=False)
        resP, _ = surP.auto_run(SAMPLING, POINTS, SURRO_POLYNOM, run_validation=VALIDATION)
        POINTS = 14
        surRg = Surrogate(use_abaqus=True, pgf=PGF, show_plots=SHOW_PLOT, scale_it=True)
        resRg, _ = surRg.auto_run(SAMPLING, POINTS, SURRO_RBF, run_validation=VALIDATION, params=[1., 'gaus'])
        POINTS = 20
        surRm = Surrogate(use_abaqus=True, pgf=PGF, show_plots=SHOW_PLOT, scale_it=True)
        resRm, _ = surRm.auto_run(SAMPLING, POINTS, SURRO_RBF, run_validation=VALIDATION, params=[1., 'multi-quadratic'])
        POINTS = 16
        surK = Surrogate(use_abaqus=True, pgf=PGF, show_plots=SHOW_PLOT, scale_it=False)
        resK, _ = surK.auto_run(SAMPLING, POINTS, SURRO_KRIGING, run_validation=VALIDATION)

        print('POLY------------------------')
        surP.print_results()
        print('RBF gaus------------------------')
        surRg.print_results()
        print('RBF mq------------------------')
        surRm.print_results()
        print('KRIGING------------------------')
        surK.print_results()

        #plot opti-lines
        opti_lines = PlotHelper(['Rippen', 'opt. Gewicht'], pgf=PGF, fancy=True)
        l0 = opti_lines.ax.plot(resP.opti_curve[0], resP.opti_curve[3], '-', label='Polynom') #, color='dodgerblue')
        opti_lines.ax.plot([resP.optimum_rib], [resP.optimum_weight], 'o', color=l0[0].get_color())

        l1 = opti_lines.ax.plot(resRg.opti_curve[0], resRg.opti_curve[3], '-', label='RBF-gauß') #, color='orange')
        opti_lines.ax.plot([resRg.optimum_rib], [resRg.optimum_weight], 'o', color=l1[0].get_color())

        l2 = opti_lines.ax.plot(resRm.opti_curve[0], resRm.opti_curve[3], '-', label='RBF-mq') #, color='magenta')
        opti_lines.ax.plot([resRm.optimum_rib], [resRm.optimum_weight], 'o', color=l2[0].get_color())

        l3 = opti_lines.ax.plot(resK.opti_curve[0], resK.opti_curve[3], '-', label='Kriging') #, color='mediumseagreen')
        opti_lines.ax.plot([resK.optimum_rib], [resK.optimum_weight], 'o', color=l3[0].get_color())

        opt = opti_lines.ax.plot([17], [405.5], 'kx', label='exakte Lösung')

        import matplotlib.ticker as ticker
        opti_lines.ax.xaxis.set_major_locator(ticker.IndexLocator(base=2, offset=0))
        opti_lines.finalize(height=2.3, legendLoc=9, bbox_to_anchor=(0.5, -0.42), legendNcol=3)
        opti_lines.ax.locator_params(nbins=4, axis='y')
        import matplotlib.pyplot as plt
        plt.subplots_adjust(bottom=0.49, top=0.98)
        opti_lines.save(Constants().PLOT_PATH + 'optiLinesHalton.pdf')
        opti_lines.show()