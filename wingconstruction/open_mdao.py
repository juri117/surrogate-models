__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :runs openMDAO optimization on wing-structure
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import sys
import os
from datetime import datetime
import numpy as np
import math

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/OpenMDAO')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyDOE2')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyoptsparse')
from openmdao.api import Problem, ExecComp, pyOptSparseDriver, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingconstruction.wingutils.constants import Constants
from wingconstruction.multi_run import MultiRun
from wingconstruction.wingutils.defines import *
from myutils.plot_helper import PlotHelper
from myutils.time_track import TimeTrack

PROJECT_NAME_PREFIX = 'iterSLSQP_FEM'

LOG_FILE_PATH = Constants().WORKING_DIR + '/' + PROJECT_NAME_PREFIX + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'

WEIGHT_PANALTY_FAC = 0
USE_ABA = True
PGF = False
USE_SCALING = True

offset_rib = range_rib[0]
offset_shell = range_shell[0]
scale_rib = range_rib[1] - range_rib[0]
scale_shell = range_shell[1] - range_shell[0]

offset_weight = 0.
offset_stress = 0.
scale_weight = 1.
scale_stress = 1.


class WingStructure(ExplicitComponent):

    def setup(self):
        ######################
        ### needed Objects ###
        self.runner = MultiRun(use_calcu=not USE_ABA, use_aba=USE_ABA, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX, force_recalc=False)

        #####################
        ### openMDAO init ###
        ### INPUTS
        self.add_input('ribs', val=int((22 - offset_rib) / scale_rib), desc='number of ribs')
        self.add_input('shell', val=(0.003 - offset_shell) / scale_shell, desc='thickness of the shell')

        ### OUTPUTS
        self.add_output('stress', val=1e8)#, ref=1e8)
        self.add_output('weight', val=100.)#, ref=100)

        self.declare_partials('*', '*', method='fd')
        self.executionCounter = 0
        self.timer = TimeTrack()

    def compute(self, inputs, outputs):
        ribs = (inputs['ribs'][0] * scale_rib) + offset_rib
        rib0 = int(math.floor(ribs))
        rib1 = int(math.ceil(ribs))

        ribs_int = int(round(ribs))

        shell = (inputs['shell'][0] * scale_shell) + offset_shell
        self.runner.project_name_prefix = PROJECT_NAME_PREFIX + '_{:05d}'.format(self.executionCounter)

        pro0 = self.runner.new_project_r_t(rib0, shell)
        pro0 = self.runner.run_project(pro0, used_cpus=1)
        res0 = pro0.resultsCalcu
        pro1 = self.runner.new_project_r_t(rib1, shell)
        pro1 = self.runner.run_project(pro1, used_cpus=1)
        res1 = pro1.resultsCalcu

        if USE_ABA:
            res0 = pro0.resultsAba
            res1 = pro1.resultsAba

        if rib1 - rib0 < 0.000000001:
            stress = res0.stressMisesMax
            weight = pro0.calc_wight()
        else:
            stress = (res0.stressMisesMax + ((ribs) - rib0)
                      * ((res1.stressMisesMax - res0.stressMisesMax) / (rib1 - rib0)))
            weight = (pro0.calc_wight() + ((ribs) - rib0)
                      * ((pro1.calc_wight() - pro0.calc_wight()) / (rib1 - rib0)))

        outputs['stress'] = (stress - offset_stress) / scale_stress
        outputs['weight'] = (weight - offset_weight) / scale_weight

        if WEIGHT_PANALTY_FAC > 0:
            weight_panalty = ((ribs) % 1)
            if weight_panalty >= 0.5:
                weight_panalty = 1. - weight_panalty
            outputs['weight'] = weight + ((weight_panalty * WEIGHT_PANALTY_FAC) * WEIGHT_FAC)

        write_mdao_log(str(self.executionCounter) + ','
                     + str(self.timer.get_time()) + ','
                     + str(ribs) + ','
                     + str(ribs_int) + ','
                     + str(shell) + ','
                     + str(stress) + ','
                     + str(weight))
        self.executionCounter += 1
        print(
            '#{:d}: {:0.10f}({:d}), {:0.10f} -> {:0.10f}, {:0.10f}'.format(self.executionCounter, ribs, ribs_int, shell,
                                                                           stress, weight))
        print('#{:d}: {:0.10f}, {:0.10f} -> {:0.10f}, {:0.10f}'.format(self.executionCounter, inputs['ribs'][0],
                                                                       inputs['shell'][0], outputs['stress'][0],
                                                                       outputs['weight'][0]))


def write_mdao_log(out_str):
    out_str = out_str.replace('[', '')
    out_str = out_str.replace(']', '')
    output_f = open(LOG_FILE_PATH, 'a') #'a' so we append the file
    output_f.write(out_str + '\n')
    output_f.close()


def run_open_mdao():
    if USE_SCALING:
        # prepare scaling
        global offset_weight
        global offset_stress
        global scale_weight
        global scale_stress
        runner = MultiRun(use_calcu=not USE_ABA, use_aba=USE_ABA, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX,
                          force_recalc=False)
        p1 = runner.new_project_r_t(range_rib[0], range_shell[0])
        runner.run_project(p1)
        offset_weight = p1.calc_wight()
        max_stress = p1.resultsAba.stressMisesMax
        p2 = runner.new_project_r_t(range_rib[1], range_shell[1])
        runner.run_project(p2)
        max_weight = p2.calc_wight()
        offset_stress = p2.resultsAba.stressMisesMax
        scale_weight = (max_weight - offset_weight)
        scale_stress = (max_stress - offset_stress)

    write_mdao_log('iter,time,ribs(float),ribs,shell,stress,weight')

    model = Group()

    indeps = IndepVarComp()

    indeps.add_output('ribs', (20 - offset_rib) / scale_rib)
    indeps.add_output('shell', (0.0025 - offset_shell) / scale_shell)

    model.add_subsystem('des_vars', indeps)
    model.add_subsystem('wing', WingStructure())
    model.connect('des_vars.ribs', ['wing.ribs', 'con_cmp1.ribs'])
    model.connect('des_vars.shell', 'wing.shell')

    # design variables, limits and constraints
    model.add_design_var('des_vars.ribs', lower=(range_rib[0] - offset_rib) / scale_rib,
                         upper=(range_rib[1] - offset_rib) / scale_rib)
    model.add_design_var('des_vars.shell', lower=(range_shell[0] - offset_shell) / scale_shell,
                         upper=(range_shell[1] - offset_shell) / scale_shell)

    # objective
    model.add_objective('wing.weight', scaler=1)

    # constraint
    print('constrain stress: ' + str((max_shear_strength - offset_stress) / scale_stress))
    model.add_constraint('wing.stress', upper=(max_shear_strength - offset_stress) / scale_stress)
    model.add_subsystem('con_cmp1',
                        ExecComp('con1 = (ribs * ' + str(scale_rib) + ') - int(ribs[0] * ' + str(scale_rib) + ')'))
    model.add_constraint('con_cmp1.con1', upper=.5)

    prob = Problem(model)

    # setup the optimization
    if USE_PYOPTSPARSE:
        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = OPTIMIZER
        prob.driver.opt_settings['SwarmSize'] = 6
        prob.driver.opt_settings['stopIters'] = 5
    else:
        prob.driver =  ScipyOptimizeDriver()
        prob.driver.options['optimizer'] = OPTIMIZER  # ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
        prob.driver.options['tol'] = TOL
        prob.driver.options['disp'] = True
        prob.driver.options['maxiter'] = 3
        #prob.driver.opt_settings['etol'] = 100

    prob.setup()
    prob.set_solver_print(level=0)
    prob.model.approx_totals()
    prob.setup(check=True, mode='fwd')
    prob.run_driver()

    print('done')
    print('ribs: ' + str((prob['wing.ribs'] * scale_rib) + offset_rib))
    print('shell: ' + str((prob['wing.shell'] * scale_shell) + offset_shell) + ' m')
    print('weight= ' + str((prob['wing.weight'] * scale_weight) + offset_weight))
    print('stress= ' + str((prob['wing.stress'] * scale_stress) + offset_stress) + ' ~ ' + str(prob['wing.stress']))
    print('execution counts wing: ' + str(prob.model.wing.executionCounter))


def plot_iter(file_path=None):
    if file_path == None:
        file_path = LOG_FILE_PATH
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    iter = data[:, 0]
    time = data[:, 1]
    ribs = data[:, 2]
    shell = data[:, 4]
    stress = data[:, 5]
    weight = data[:, 6]
    #print some info:
    print('number of iterations: {:d}'.format(int(iter[-1])+1))
    print('total time: {:f}'.format(time[-1] - time[0]))
    iter_plot = PlotHelper([], fancy=True, pgf=PGF)
    ax1 = iter_plot.fig.add_subplot(211)
    ax2 = iter_plot.fig.add_subplot(212)
    # param plot
    iter_param = PlotHelper(['', 'Rippen'], fancy=True, ax=ax1, pgf=PGF)
    iter_param.ax.plot(iter, ribs, color='teal')
    iter_param.ax.yaxis.label.set_color('teal')
    ax_shell = iter_param.ax.twinx()
    ax_shell.set_ylabel('Blechd. in mm')
    ax_shell.yaxis.label.set_color('orange')
    ax_shell.plot(iter, shell * 1000, color='orange')
    iter_param.ax.set_ylim(range_rib)
    ax_shell.set_ylim(tuple(1000*x for x in range_shell))
    iter_param.finalize(show_legend=False)
    from matplotlib.ticker import MaxNLocator
    iter_param.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    # results plot
    iter_res = PlotHelper(['Iteration', 'Mises in Pa'], fancy=True, ax=ax2, pgf=PGF)
    iter_res.ax.plot(iter, stress, color='tomato')
    iter_res.ax.plot([min(iter), max(iter)], [max_shear_strength, max_shear_strength], '--', color='tomato', label='max. Spannung')
    iter_res.ax.yaxis.label.set_color('tomato')
    ax_weight = iter_res.ax.twinx()
    ax_weight.set_ylabel('Gewicht in kg')
    ax_weight.yaxis.label.set_color('royalblue')
    ax_weight.plot(iter, weight, color='royalblue')
    iter_res.ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    #leg = iter_res.finalize(legendLoc='upper right', show_legend=True, bbox_to_anchor=(.95, 1.25))
    iter_plot.finalize(width=6, height=3, tighten_layout=True)
    handles, labels = iter_res.ax.get_legend_handles_labels()
    leg = iter_plot.fig.legend(handles, labels, loc='lower center', bbox_to_anchor=(0.65, 0.42), ncol=2, fancybox=True)
    leg.get_frame().set_facecolor('#FFFFFF')
    leg.get_frame().set_alpha(1.)
    import matplotlib.pyplot as plt
    plt.subplots_adjust(wspace=0.025, hspace=0.55, bottom=0.18, top=.98)
    #iter_param.ax.text(7.5, 25, 'Eingägne')
    #iter_res.ax.text(7.5, 6.5e+8, 'Ausgägne')
    iter_plot.save('../data_out/plot/openMDAOconv_ALPSO.pdf', transparent=False)
    iter_plot.show()


if __name__ == '__main__':
    #01
    SHELL_FACTOR = 1e-2  # 1e-2
    RIB_FACTOR = 1e-6  # 1e-6
    WEIGHT_FAC = 1e-3
    STRESS_FAC = 1e-8
    TOL = 1e-2
    USE_PYOPTSPARSE = False
    OPTIMIZER = 'SLSQP'
    WEIGHT_PANALTY_FAC = 0

    #02
    if False:
        SHELL_FACTOR = 1  # 1e-2
        RIB_FACTOR = 1e-6  # 1e-6
        WEIGHT_FAC = 1e-3
        STRESS_FAC = 1e-8
        TOL = 1e-3
        USE_PYOPTSPARSE = True
        OPTIMIZER = 'ALPSO'
        WEIGHT_PANALTY_FAC = 10


    run_open_mdao()
    plot_iter()

    #plot_iter(Constants().WORKING_DIR + '/' + 'iterSLSQP_FEM2018-10-10_16_10_59_V001.csv')
    #plot_iter(Constants().WORKING_DIR + '/' + 'iterALPSO_FEM2018-10-21_11_00_27_V002.csv')