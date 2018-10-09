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
from openmdao.api import pyOptSparseDriver, Problem, ExecComp, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingConstruction.wingUtils.Constants import Constants
from wingConstruction.MultiRun import MultiRun
from wingConstruction.wingUtils.defines import *
from myUtils.PlotHelper import PlotHelper

PROJECT_NAME_PREFIX = 'v2iterSLSQP'

LOG_FILE_PATH = Constants().WORKING_DIR + '/' + PROJECT_NAME_PREFIX + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'

USE_ABA = True

class WingStructureV2(ExplicitComponent):

    def setup(self):
        ######################
        ### needed Objects ###
        self.runner = MultiRun(use_calcu=not USE_ABA, use_aba=USE_ABA, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX, force_recalc=False)

        #####################
        ### openMDAO init ###
        ### INPUTS
        self.add_input('ribs', val=int((range_rib[0] + range_rib[1]) / 2), desc='number of ribs')
        self.add_input('shell', val=((range_shell[0] + range_shell[1]) / 2), desc='thickness of the shell')

        ### OUTPUTS
        self.add_output('stress', val=1e8)#, ref=1e8)
        self.add_output('weight', val=100.)#, ref=100)

        self.declare_partials('*', '*', method='exact')
        self.executionCounter = 0

    def compute(self, inputs, outputs):
        rib0 = int(math.floor(inputs['ribs'][0]))
        rib1 = int(math.ceil(inputs['ribs'][0]))
        ribs = int(round(inputs['ribs'][0]))
        shell = inputs['shell'][0]
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
            stress = (res0.stressMisesMax + ((inputs['ribs'][0]) - rib0)
                      * ((res1.stressMisesMax - res0.stressMisesMax) / (rib1 - rib0)))
            weight = (pro0.calc_wight() + ((inputs['ribs'][0]) - rib0)
                      * ((pro1.calc_wight() - pro0.calc_wight()) / (rib1 - rib0)))

        outputs['stress'] = stress
        outputs['weight'] = weight

        write_to_log(str(self.executionCounter) + ','
                     + datetime.now().strftime('%H:%M:%S') + ','
                     + str(inputs['ribs']) + ','
                     + str(ribs) + ','
                     + str(inputs['shell']) + ','
                     + str(outputs['stress']) + ','
                     + str(outputs['weight']))
        self.executionCounter += 1
        print('#{:d}: {:0.10f}({:d}), {:0.10f} -> {:0.10f}, {:0.10f}'.format(self.executionCounter, inputs['ribs'][0], ribs, inputs['shell'][0], outputs['stress'][0], outputs['weight'][0]))


def write_to_log(out_str):
    out_str = out_str.replace('[', '')
    out_str = out_str.replace(']', '')
    output_f = open(LOG_FILE_PATH, 'a') #'a' so we append the file
    output_f.write(out_str + '\n')
    output_f.close()


def run_open_mdao():
    prob = Problem()
    indeps = prob.model.add_subsystem('indeps', IndepVarComp())
    indeps.add_output('ribs', 14.)
    indeps.add_output('shell', 0.02)
    prob.model.add_subsystem('wing', WingStructureV2())
    prob.model.connect('indeps.ribs', 'wing.ribs')
    prob.model.connect('indeps.shell', 'wing.shell')

    '''
    ALPSO - Augmented Lagrangian Particle Swarm Optimizer
    
    '''

    prob.driver = pyOptSparseDriver()
    prob.driver.options['optimizer'] = 'ALPSO'
    #prob.driver.options['tol'] = 1e-6
    #prob.driver.options['maxiter'] = 100000
    #prob.driver.options['tol'] = 1e-9
    #prob.driver.options['disp'] = True
    prob.model.add_design_var('indeps.ribs', lower=5, upper=30)
    prob.model.add_design_var('indeps.shell', lower=0.006, upper=0.03)
    prob.model.add_objective('wing.weight', scaler=100)
    prob.model.add_constraint('wing.stress', upper=5.72e8 / 1.5, scaler=1e8)

    prob.setup()
    prob.run_driver()

    print('done')
    print('ribs: ' + str(prob['wing.ribs']))
    print('shell: ' + str(prob['wing.shell']) + ' m')
    print('weight= ' + str(prob['wing.weight']))
    print('stress= ' + str(prob['wing.stress']))
    print('execution counts wing: ' + str(prob.model.wing.executionCounter))


def plot_iter(file_path=None):
    if file_path == None:
        file_path = LOG_FILE_PATH
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    iter = data[:, 0]
    ribs = data[:, 2]
    shell = data[:, 4]
    stress = data[:, 5]
    weight = data[:, 6]
    iter_plot = PlotHelper([], fancy=False, pgf=False)
    ax1 = iter_plot.fig.add_subplot(211)
    ax2 = iter_plot.fig.add_subplot(212)
    # param plot
    iter_param = PlotHelper(['', 'Rippen'], fancy=False, ax=ax1, pgf=False)
    iter_param.ax.plot(iter, ribs, color='teal')
    iter_param.ax.yaxis.label.set_color('teal')
    ax_shell = iter_param.ax.twinx()
    ax_shell.set_ylabel('Shell')
    ax_shell.yaxis.label.set_color('orange')
    ax_shell.plot(iter, shell, color='orange')
    iter_param.finalize(width=6, height=2.5, show_legend=False)
    # results plot
    iter_res = PlotHelper(['Iteration', 'Mises in Pa'], fancy=False, ax=ax2, pgf=False)
    iter_res.ax.plot(iter, stress, color='tomato')
    iter_res.ax.plot([min(iter), max(iter)], [max_shear_strength, max_shear_strength], '--', color='tomato', label='max. Spannung')
    iter_res.ax.yaxis.label.set_color('tomato')
    ax_weight = iter_res.ax.twinx()
    ax_weight.set_ylabel('Gewicht in kg')
    ax_weight.yaxis.label.set_color('royalblue')
    ax_weight.plot(iter, weight, color='royalblue')
    iter_res.finalize(width=6, height=2.5, legendLoc='lower center', show_legend=True)

    iter_plot.save('../dataOut/plot/openMDAOconv.pdf')
    iter_plot.show()

if __name__ == '__main__':
    run_open_mdao()
    plot_iter()
    #plot_iter(file_path=Constants().WORKING_DIR + '/' + 'iterSLSQP_2018-10-03_10_48_20_aba.csv')