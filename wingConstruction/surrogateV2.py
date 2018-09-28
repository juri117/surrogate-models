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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyKriging')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/inspyred')
from pyKriging.krige import kriging as PyKriging

import numpy as np
from scipy import optimize

from wingConstruction.wingUtils.Constants import Constants
from wingConstruction.MultiRun import MultiRun
from myUtils.TimeTrack import TimeTrack
from myUtils.PlotHelper import PlotHelper
from myLibs.Kriging import Kriging
from myLibs.RBF import RBF
from myLibs.Polynomial import Polynomial
from myLibs.LatinHyperCube import LatinHyperCube
from myLibs.Halton import Halton
from myLibs.interface.OptiLatinHyper import OptiLatinHyper
from myLibs.StructuredSample import StructuredSample
from myLibs.Validation import Validation
from myLibs.Validation import ValidationResults
from wingConstruction.fem.WingConstructionV4 import WingConstruction
from wingConstruction.wingUtils.defines import *



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
def surrogate_analysis(sampling_type, sample_point_count, surro_type, use_abaqus=False, pgf=False, show_plots=True, force_recalc=False):
    timer = TimeTrack('surrogateAnalysis')
    timer.tic()
    results = SurroResults()
    RESULTS_FILE = '/2drun_2018-08-23_16_49_18_final01_cruiseLoad.csv'

    multi = MultiRun(use_calcu=True, use_aba=True, non_liner=False, force_recalc=force_recalc)

    ##################################################
    # collect data

    ribs, shell, stress, disp, weight = multi.read_data_file(RESULTS_FILE, use_abaqus=use_abaqus)

    # input 0 is ribs
    # input 1 is shell_thick

    # reduce data range
    i_rib_0 = 0
    i_rib_n = len(ribs)-1
    for i in range(0, len(ribs)):
        if range_rib[0] >= ribs[i]:
            i_rib_0 = i
        if range_rib[1] <= ribs[i]:
            i_rib_n = i
            break
    i_shell_0 = 0
    i_shell_n = len(shell)-1
    for i in range(0, len(shell)):
        if range_shell[0] >= shell[i]:
            i_shell_0 = i
        if range_shell[1] <= shell[i]:
            i_shell_n = i
            break

    new_rib = []
    new_shell = []
    new_stress = np.zeros((i_shell_n-i_shell_0+1, i_rib_n-i_rib_0+1))
    #newWeight = np.zeros((i_shell_n-i_shell_0, i_rib_n-i_rib_0))
    for r in range(0, i_rib_n-i_rib_0+1):
        new_rib.append(ribs[i_rib_0+r])
        for s in range(0, i_shell_n-i_shell_0+1):
            new_stress[s][r] = stress[i_shell_0+s][i_rib_0+r]
            #newWeight[s][r] = weight[i_shell_0 + s][i_rib_0 + r]

    for s in range(0, i_shell_n-i_shell_0+1):
        new_shell.append(shell[i_shell_0+s])

    ribs = new_rib
    shell = new_shell
    stress = new_stress
    #weight = newWeight

    ##################################################
    # sample plan

    if sampling_type == SAMPLE_LATIN:
        sam = LatinHyperCube()
    elif sampling_type == SAMPLE_HALTON:
        sam = Halton()
    elif sampling_type == SAMPLE_STRUCTURE:
        sam = StructuredSample()
    elif sampling_type == SAMPLE_OPTI_LATIN_HYPER:
        sam = OptiLatinHyper()
    else:
        print('unknown sample plan selected')
        results.errorStr = 'unknown sample plan selected'
        return results

    sample_points = sam.generate_sample_plan(sample_point_count, 2, [range_rib, range_shell])
    # make the ribs be int
    for i in range(0, len(sample_points)):
        sample_points[i][0] = int(round(sample_points[i][0]))
    known_params = np.array(sample_points)
    known_rib = known_params[:, 0]
    known_shell = known_params[:, 1]

    print('sample plan using {:d} known values'.format(len(known_rib)))

    ##################################################
    # FEM calculation, collecting results

    known_stress = multi.run_sample_points(known_rib, known_shell, use_abaqus=use_abaqus)

    ##################################################
    # build surrogate model and fit it

    # will be needed later for Validation
    update_params = None

    if surro_type == SURRO_KRIGING:
        surro_class = Kriging
        surro = Kriging(known_params, known_stress)
        # prev stored results:
        surro.update_param([0.002261264770141511, 277826.21903867245], [1.8766170168043503, 1.9959876593551822])
        print('starting Likelihood optimization')
        surro.optimize()
        if show_plots:
            pltLike = surro.plot_likelihoods(pgf=pgf)
            pltLike.save(Constants().PLOT_PATH + 'wingSurroLikely.pdf')
            minLike = surro.calc_likelihood()
            print('minLike = ' + str(minLike))
            print('@theta1 = ' + str(surro.get_theta()[0]))
            print('@theta2 = ' + str(surro.get_theta()[1]))
            print('@p1 = ' + str(surro.get_p()[0]))
            print('@p2 = ' + str(surro.get_p()[1]))
    elif surro_type == SURRO_RBF:
        surro_class = RBF
        surro = RBF(known_params, known_stress)
        a = -.06 #.078
        surro.update_param(a, 'multi-quadratic')
        update_params = [a, 'multi-quadratic']
        if show_plots:
            print('coeff1 = ' + str(surro.get_coeff()[0]))
            print('coeff2 = ' + str(surro.get_coeff()[1]))
    elif surro_type == SURRO_POLYNOM:
        surro_class = Polynomial
        surro = Polynomial(known_params, known_stress)
        o = 3
        surro.update_param(o)
        surro.generate_formula()
        update_params = [o]
    elif surro_type == SURRO_PYKRIGING:
        surro_class = PyKriging
        surro = PyKriging(known_params, known_stress)
        surro.train()
    else:
        print('unknown surrogate type selected')
        results.errorStr = 'unknown surrogate type selected'
        return results

    ##################################################
    # optimize

    def shell_predict(shell_thick, surro_inst, rib_num):
        #krig_inst = args[0]
        #rib_num = args[1]
        stress_val = surro_inst.predict([rib_num, shell_thick])
        return stress_val - max_shear_strength

    opti_ribs = []
    opti_shell = []
    opti_stress = []
    opti_weights = []
    used_ribs = list(range(range_rib[0], range_rib[1]+1))#list(range(int(min(known_rib)), int(max(known_rib))))
    for i in range(0, len(used_ribs)):
        # SLSQP: proplem; find local min not glob. depending on init-vals
        init_guess = (min(known_shell) + max(known_shell))/2
        bnds = [(min(known_shell), max(known_shell))]
        #res = minimize(shell_predict, init_guess, args=[krig, ribs[i]], method='SLSQP', tol=1e-6, options={'disp': True, 'maxiter': 99999}, bounds=bnds)
        #opti_shell.append(res.x[0])
        root = optimize.newton(shell_predict, init_guess, args=[surro, used_ribs[i]])
        opti_ribs.append(used_ribs[i])
        opti_shell.append(root)
        opti_stress.append(surro.predict([used_ribs[i], root]))
        weight = WingConstruction.calc_weight_stat(wing_length,
                                                   chord_length * 0.4,
                                                   chord_height,
                                                   used_ribs[i],
                                                   root,
                                                   density)
        opti_weights.append(weight)
        """
        plt.clf()
        sll = np.linspace(min(shell), max(shell), 500)
        strs = []
        for s in sll:
            strs.append(shell_predict(s, krig, ribs[i]))
        plt.plot(sll, strs, 'b-')
        opt_strs = shell_predict(root, krig, ribs[i])
        plt.plot([root], [opt_strs], 'ro')
        plt.grid(True)
        plt.show()
        """

    # exclude model edges from opti vals
    opti_ribs = opti_ribs[:-1]
    opti_ribs = opti_ribs[1:]
    opti_shell = opti_shell[:-1]
    opti_shell = opti_shell[1:]
    opti_stress = opti_stress[:-1]
    opti_stress = opti_stress[1:]
    opti_weights = opti_weights[:-1]
    opti_weights = opti_weights[1:]

    best_i = opti_weights.index(min(opti_weights))

    if show_plots:
        optWeightPlot = PlotHelper(['ribs', 'weight'], pgf=pgf)
        optWeightPlot.ax.plot(opti_ribs, opti_weights, 'b-')
        optWeightPlot.ax.plot([opti_ribs[best_i]], opti_weights[best_i], 'rx', label='minimum')
        optWeightPlot.finalize()

    print('optimum:')
    print('ribs: {:f}'.format(opti_ribs[best_i]))
    print('shell: {:f}'.format(opti_shell[best_i]))
    print('weight: {:f}'.format(opti_weights[best_i]))
    print('stress: {:f}'.format(max_shear_strength))

    results.optimumRib = opti_ribs[best_i]
    results.optimumShell = opti_shell[best_i]
    results.optimumWeights = opti_weights[best_i]

    ##################################################
    # validate

    # validate points 0.002, 0.0033
    valiParams = np.array([[7., 0.0025], [13., 0.0030], [16., 0.0028], [results.optimumRib, results.optimumShell]])
    valiValues = multi.run_sample_points(valiParams.T[0], valiParams.T[1], use_abaqus=use_abaqus)
    # valiValues = np.array(list(map(f_2D, valiParams)))
    # valiParams = valiParams.reshape((len(valiParams), 1))

    p_x, p_y = np.meshgrid(ribs, shell)
    params = np.array([p_x.flatten(), p_y.flatten()]).T
    values = stress.flatten()

    vali = Validation()
    vali_r = vali.run_full_analysis(params, values,
                                    known_params, known_stress,
                                    valiParams, valiValues,
                                    surro.predict, surro_class, update_params=update_params)
    results.valiResults = vali_r
    print('avg deviation: {:.3e} (-> {:.3f}%)'.format(vali_r.deviation, vali_r.deviation * 100.))
    print('rmse: {:f}'.format(vali_r.rmse))
    print('mae: {:f}'.format(vali_r.mae))
    print('rae: {:s}'.format(str(vali_r.rae)))
    print('press: {:f}'.format(vali_r.press))

    ##################################################
    # plot it

    if show_plots:
        # convert all to np arrays
        shell = np.array(shell)
        opti_shell = np.array(opti_shell)
        known_shell = np.array(known_shell)

        plot3d = PlotHelper(['ribs', 'shell thickness in mm', 'mises stress'], fancy=False, font_size=16, pgf=pgf)

        #realDat = plot3d.ax.plot_wireframe(rib_mat, shell_mat, stress, color='g', alpha=0.5, label='fem data')

        # plot FEM data as lines
        # plot FEM data as lines

        for i in range(0,len(ribs)):
            if i == 0:
                plot3d.ax.plot(np.ones((len(shell)))*ribs[i], shell*1000., stress[:,i], 'g-', lw=3., label='fem data')
            else:
                plot3d.ax.plot(np.ones((len(shell))) * ribs[i], shell*1000., stress[:, i], 'g-', lw=3.)


        #realDatMark = plot3d.ax.scatter(rib_mat, shell_mat, stress, c='g', marker='x', label='fem measurements')

        # plot limit load as wireframe
        #limit = np.full((n_thick, n_rib),max_shear_strength)
        #plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.2, label='limit load')

        # plot surrogate model as wireframe
        ribs_sample = np.linspace(range_rib[0], range_rib[1], 200)
        shell_sample = np.linspace(range_shell[0], range_shell[1], 200)
        krigPlot = plot3d.plot_function_3d(surro.predict, ribs_sample, shell_sample, r'$\widehat{f}_{krig}$', color='b', scale=[1., 1000., 1.])
        samplePoints = plot3d.ax.plot(known_rib, known_shell*1000., known_stress, 'bo', label='sampling points')

        # plot limit load line
        plot3d.ax.plot(opti_ribs, opti_shell*1000., opti_stress, 'k--', lw=3., label='max. stress line')
        # plot optimal point
        plot3d.ax.plot([opti_ribs[best_i]], [opti_shell[best_i]*1000.], [opti_stress[best_i]], 'rx', markersize=12, markeredgewidth=5, label='global optimum')
        plot3d.ax.locator_params(nbins=7, axis='y')

        plot3d.ax.set_zlim3d(np.min(np.array(stress)), max_shear_strength*1.2)
        plot3d.ax.set_ylim3d(np.min(np.array(shell))*1000.,np.max(np.array(shell))*1000.)

        plot3d.finalize(height=7, width=9, legendLoc=8, legendNcol=3, bbox_to_anchor=(0.5, -0.25), tighten_layout=True)
        plot3d.ax.view_init(18, 40)
        plot3d.save(Constants().PLOT_PATH + 'wingSurro.pdf')
        plot3d.show()
    results.runtime = timer.toc()
    return results


class SurroResults:

    def __init__(self):
        self.optimumRib = 0.
        self.optimumShell = 0.
        self.optimumWeights = 0.
        self.runtime = 0.
        self.errorStr = '-'
        self.valiResults = ValidationResults()
        #self.deviation = 0.
        #self.rmse = 0.
        #self.mae = 0.
        #self.rae = 0.
        #self.press = 0.


if __name__ == '__main__':
    # SAMPLE_LATIN, SAMPLE_HALTON, SAMPLE_STRUCTURE, SAMPLE_OPTI_LATIN_HYPER
    # SURRO_KRIGING, SURRO_RBF, SURRO_POLYNOM, SURRO_PYKRIGING
    surrogate_analysis(SAMPLE_LATIN, 14, SURRO_POLYNOM, use_abaqus=False, pgf=False, show_plots=True)