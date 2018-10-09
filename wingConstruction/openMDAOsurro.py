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
from openmdao.api import Problem, pyOptSparseDriver, ExecComp, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingConstruction.wingUtils.Constants import Constants
from wingConstruction.MultiRun import MultiRun
from wingConstruction.wingUtils.defines import *
from myUtils.PlotHelper import PlotHelper
from wingConstruction.surrogateV2 import surrogate_analysis
from myUtils.TimeTrack import TimeTrack

PROJECT_NAME_PREFIX = 'iterSLSQP'

LOG_FILE_PATH = Constants().WORKING_DIR + '/' + PROJECT_NAME_PREFIX + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
SHELL_FACTOR = 1#1e-2
RIB_FACTOR = 1#1e-6
WEIGHT_FAC = 1e-3
STRESS_FAC = 1e-8

#WEIGHT_PANALTY_FAC = 10.

USE_ABA = True

PGF = True

class WingStructureSurro(ExplicitComponent):

    def setup(self):
        ######################
        ### needed Objects ###
        self.runner = MultiRun(use_calcu=not USE_ABA, use_aba=USE_ABA, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX, force_recalc=False)
        _, self.surro = surrogate_analysis(SAMPLE_LATIN, 14, SURRO_POLYNOM, use_abaqus=USE_ABA, pgf=False, show_plots=False)


        #####################
        ### openMDAO init ###
        ### INPUTS
        self.add_input('ribs', val=int((range_rib[0] + range_rib[1]) / 2)*RIB_FACTOR, desc='number of ribs')
        self.add_input('shell', val=((range_shell[0] + range_shell[1]) / 2)*SHELL_FACTOR, desc='thickness of the shell')

        ### OUTPUTS
        self.add_output('stress', val=1e8)#, ref=1e8)
        self.add_output('weight', val=100.)#, ref=100)

        self.declare_partials('*', '*', method='fd')
        self.executionCounter = 0
        self.timer = TimeTrack()

    def compute(self, inputs, outputs):
        rib0 = int(math.floor(inputs['ribs'][0] / RIB_FACTOR))
        rib1 = int(math.ceil(inputs['ribs'][0] / RIB_FACTOR))

        ribs = round(inputs['ribs'][0] / RIB_FACTOR)

        shell = inputs['shell'][0] / SHELL_FACTOR
        #self.runner.project_name_prefix = PROJECT_NAME_PREFIX + '_{:05d}'.format(self.executionCounter)

        #pro = self.runner.new_project_r_t(ribs, shell)
        #pro = self.runner.run_project(pro)
        #stress = pro.resultsCalcu.stressMisesMax * STRESS_FAC
        #weight = pro.calc_wight() * WEIGHT_FAC

        pro0 = self.runner.new_project_r_t(rib0, shell)
        pro1 = self.runner.new_project_r_t(rib1, shell)

        if rib1 - rib0 < 0.000000001:
            weight = pro0.calc_wight() * WEIGHT_FAC
        else:
            weight = (pro0.calc_wight() + ((inputs['ribs'][0] / RIB_FACTOR) - rib0)
                  * ((pro1.calc_wight() - pro0.calc_wight()) / (rib1 - rib0))) \
                  * WEIGHT_FAC

        stress = self.surro.predict([ribs, shell]) * STRESS_FAC

        #weight_penalty = 0.
        #if stress > max_shear_strength * STRESS_FAC:
        #    weight_penalty = (stress - (max_shear_strength * STRESS_FAC)) * 100

        outputs['stress'] = stress
        outputs['weight'] = weight #+ (weight_penalty * WEIGHT_FAC)

        #weight_panalty = ((inputs['ribs'][0] / RIB_FACTOR) % 1)
        #if weight_panalty >= 0.5:
        #    weight_panalty = 1. - weight_panalty
        #weight_panalty = 0.

        #outputs['weight'] = (pro.calc_wight() + (weight_panalty * WEIGHT_PANALTY_FAC)) * WEIGHT_FAC

        write_mdao_log(str(self.executionCounter) + ','
                     + str(self.timer.get_time()) + ','
                     + str(inputs['ribs'] / RIB_FACTOR) + ','
                     + str(int(round(ribs))) + ','
                     + str(inputs['shell'] / SHELL_FACTOR) + ','
                     + str(outputs['stress'] / STRESS_FAC) + ','
                     + str(outputs['weight'] / WEIGHT_FAC))
        self.executionCounter += 1
        print('#{:d}: {:0.10f}({:d}), {:0.10f} -> {:0.10f}, {:0.10f}'.format(self.executionCounter, inputs['ribs'][0], int(round(ribs)), inputs['shell'][0], outputs['stress'][0], outputs['weight'][0]))




def write_mdao_log(out_str):
    out_str = out_str.replace('[', '')
    out_str = out_str.replace(']', '')
    output_f = open(LOG_FILE_PATH, 'a') #'a' so we append the file
    output_f.write(out_str + '\n')
    output_f.close()


def run_open_mdao():
    write_mdao_log('iter,time,ribs(float),ribs,shell,stress,weight')

    model = Group()

    #indeps = prob.model.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
    indeps = IndepVarComp()
    #indeps.add_output('ribs', int((range_rib[0] + range_rib[1]) / 2) * RIB_FACTOR)
    #indeps.add_output('shell', ((range_shell[0] + range_shell[1]) / 2)*SHELL_FACTOR)

    indeps.add_output('ribs', 20 * RIB_FACTOR)
    indeps.add_output('shell', 0.0025 * SHELL_FACTOR)

    model.add_subsystem('des_vars', indeps)
    model.add_subsystem('wing', WingStructureSurro())
    model.connect('des_vars.ribs', ['wing.ribs', 'con_cmp1.ribs'])
    model.connect('des_vars.shell', 'wing.shell')

    # design variables, limits and constraints
    model.add_design_var('des_vars.ribs', lower=range_rib[0]*RIB_FACTOR, upper=range_rib[1]*RIB_FACTOR)
    model.add_design_var('des_vars.shell', lower=range_shell[0]*SHELL_FACTOR, upper=range_shell[1]*SHELL_FACTOR)

    # objective
    model.add_objective('wing.weight', scaler=1)

    # constraint
    model.add_constraint('wing.stress', upper=max_shear_strength * STRESS_FAC)
    model.add_subsystem('con_cmp1', ExecComp('con1 = (ribs * '+str(RIB_FACTOR)+') - int(ribs[0] * '+str(RIB_FACTOR)+')'))
    model.add_constraint('con_cmp1.con1', upper=0.1)

    #model.add_constraint('des_vars.ribs', lower=8*RIB_FACTOR, upper=30*RIB_FACTOR)
    #model.add_constraint('des_vars.shell', lower=0.006*SHELL_FACTOR, upper=0.03*SHELL_FACTOR)

    prob = Problem(model)

    # setup the optimization
    if USE_PYOPTSPARSE:
        prob.driver = pyOptSparseDriver()
        prob.driver.options['optimizer'] = OPTIMIZER
        prob.driver.opt_settings['SwarmSize'] = 6
    else:
        prob.driver =  ScipyOptimizeDriver()
        prob.driver.options['optimizer'] = OPTIMIZER  # ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
        prob.driver.options['tol'] = TOL
        prob.driver.options['disp'] = True
        prob.driver.options['maxiter'] = 100
        prob.driver.opt_settings['etol'] = 100

    prob.setup()
    prob.set_solver_print(level=0)
    prob.model.approx_totals()
    prob.setup(check=True, mode='fwd')
    prob.run_driver()

    print('done')
    print('ribs: ' + str(prob['wing.ribs'] / RIB_FACTOR))
    print('shell: ' + str(prob['wing.shell'] / SHELL_FACTOR) + ' m')

    print('weight= ' + str(prob['wing.weight'] / WEIGHT_FAC))
    print('stress= ' + str(prob['wing.stress'] / STRESS_FAC))

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
    print('number of iterations: {:d}'.format(int(iter[-1])))
    print('total time: {:f}'.format(time[-1] - time[0]))
    iter_plot = PlotHelper([], fancy=True, pgf=PGF)
    ax1 = iter_plot.fig.add_subplot(211)
    ax2 = iter_plot.fig.add_subplot(212)
    # param plot
    iter_param = PlotHelper(['', 'Rippen'], fancy=True, ax=ax1, pgf=PGF)
    iter_param.ax.plot(iter, ribs, color='teal')
    iter_param.ax.yaxis.label.set_color('teal')
    ax_shell = iter_param.ax.twinx()
    ax_shell.set_ylabel('Shell')
    ax_shell.yaxis.label.set_color('orange')
    ax_shell.plot(iter, shell, color='orange')
    iter_param.ax.set_ylim(range_rib)
    ax_shell.set_ylim(range_shell)
    iter_param.finalize(width=6, height=2.5, show_legend=False)
    # results plot
    iter_res = PlotHelper(['Iteration', 'Mises in Pa'], fancy=True, ax=ax2, pgf=PGF)
    iter_res.ax.plot(iter, stress, color='tomato')
    iter_res.ax.plot([min(iter), max(iter)], [max_shear_strength, max_shear_strength], '--', color='tomato', label='max. Spannung')
    iter_res.ax.yaxis.label.set_color('tomato')
    ax_weight = iter_res.ax.twinx()
    ax_weight.set_ylabel('Gewicht in kg')
    ax_weight.yaxis.label.set_color('royalblue')
    ax_weight.plot(iter, weight, color='royalblue')
    iter_res.finalize(width=6, height=2.5, legendLoc='upper right', show_legend=True)

    iter_plot.save('../dataOut/plot/openMDAOconv_ALPSO.pdf')
    iter_plot.show()

if __name__ == '__main__':
    #01
    SHELL_FACTOR = 1  # 1e-2
    RIB_FACTOR = 1e-6  # 1e-6
    WEIGHT_FAC = 1e-3
    STRESS_FAC = 1e-8
    TOL = 1e-3
    USE_PYOPTSPARSE = False
    OPTIMIZER = 'SLSQP'

    #02
    if True:
        SHELL_FACTOR = 1  # 1e-2
        RIB_FACTOR = 1e-6  # 1e-6
        WEIGHT_FAC = 1e-3
        STRESS_FAC = 1e-8
        TOL = 1e-3
        USE_PYOPTSPARSE = True
        OPTIMIZER = 'ALPSO'

    run_open_mdao()
    plot_iter()