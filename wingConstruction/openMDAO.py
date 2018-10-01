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
from openmdao.api import Problem, ExecComp, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingConstruction.wingUtils.Constants import Constants
from wingConstruction.MultiRun import MultiRun
from wingConstruction.wingUtils.defines import *

LOG_FILE_PATH = Constants().WORKING_DIR + '/om_iterations_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'

PROJECT_NAME_PREFIX = 'iter'

SHELL_FACTOR = 1 #1e-2
RIB_FACTOR = 1e-6
WEIGHT_FAC = 1e-4
STRESS_FAC = 1e-9

WEIGHT_PANALTY_FAC = 10.

class WingStructure(ExplicitComponent):

    def setup(self):
        ######################
        ### needed Objects ###
        self.runner = MultiRun(use_calcu=True, use_aba=False, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX, force_recalc=False)

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

    def compute(self, inputs, outputs):
        ribTop = int(math.ceil(inputs['ribs'][0] / RIB_FACTOR))
        ribBot = int(math.floor(inputs['ribs'][0] / RIB_FACTOR))
        ribs = int(round(inputs['ribs'][0] / RIB_FACTOR))

        shell = inputs['shell'][0] / SHELL_FACTOR
        self.runner.project_name_prefix = PROJECT_NAME_PREFIX + '_{:05d}'.format(self.executionCounter)

        #pro = self.runner.new_project_r_t(ribs, shell)
        #pro = self.runner.run_project(pro)
        #stress = pro.resultsCalcu.stressMisesMax * STRESS_FAC
        #weight = pro.calc_wight() * WEIGHT_FAC

        proT = self.runner.new_project_r_t(ribTop, shell)
        proT = self.runner.run_project(proT)
        proB = self.runner.new_project_r_t(ribBot, shell)
        proB = self.runner.run_project(proB)
        if ribTop - ribBot < 0.000000001:
            stress = proT.resultsCalcu.stressMisesMax * STRESS_FAC
            weight = proT.calc_wight() * WEIGHT_FAC
        else:
            stress = (proB.resultsCalcu.stressMisesMax + (inputs['ribs'][0] - ribBot)
                        * ((proT.resultsCalcu.stressMisesMax - proB.resultsCalcu.stressMisesMax)/(ribTop - ribBot)))\
                        * STRESS_FAC
            weight = (proB.calc_wight() + (inputs['ribs'][0] - ribBot)
                        * ((proT.calc_wight() - proB.calc_wight()) / (ribTop - ribBot)))\
                        * WEIGHT_FAC

        outputs['stress'] = stress
        outputs['weight'] = weight

        #weight_panalty = ((inputs['ribs'][0] / RIB_FACTOR) % 1)
        #if weight_panalty >= 0.5:
        #    weight_panalty = 1. - weight_panalty
        #weight_panalty = 0.

        #outputs['weight'] = (pro.calc_wight() + (weight_panalty * WEIGHT_PANALTY_FAC)) * WEIGHT_FAC

        write_to_log(str(self.executionCounter) + ','
                     + datetime.now().strftime('%H:%M:%S') + ','
                     + str(inputs['ribs'] / RIB_FACTOR) + ','
                     + str(ribs) + ','
                     + str(inputs['shell'] / SHELL_FACTOR) + ','
                     + str(outputs['stress'] / STRESS_FAC) + ','
                     + str(outputs['weight'] / WEIGHT_FAC))
        self.executionCounter += 1
        print('#{:d}: {:f}({:d}), {:f} -> {:f}, {:f}'.format(self.executionCounter, inputs['ribs'][0], ribs, inputs['shell'][0], outputs['stress'][0], outputs['weight'][0]))


def write_to_log(out_str):
    out_str = out_str.replace('[', '')
    out_str = out_str.replace(']', '')
    output_f = open(LOG_FILE_PATH, 'a') #'a' so we append the file
    output_f.write(out_str + '\n')
    output_f.close()


def run_open_mdao():

    model = Group()

    #indeps = prob.model.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
    indeps = IndepVarComp()
    #indeps.add_output('ribs', int((range_rib[0] + range_rib[1]) / 2) * RIB_FACTOR)
    #indeps.add_output('shell', ((range_shell[0] + range_shell[1]) / 2)*SHELL_FACTOR)

    indeps.add_output('ribs', 10 * RIB_FACTOR)
    indeps.add_output('shell', 0.002562 * SHELL_FACTOR)

    model.add_subsystem('des_vars', indeps)
    model.add_subsystem('wing', WingStructure())
    model.connect('des_vars.ribs', ['wing.ribs', 'con_cmp1.ribs'])
    model.connect('des_vars.shell', 'wing.shell')

    # design variables, limits and constraints
    model.add_design_var('des_vars.ribs', lower=range_rib[0]*RIB_FACTOR, upper=range_rib[1]*RIB_FACTOR)
    model.add_design_var('des_vars.shell', lower=range_shell[0]*SHELL_FACTOR, upper=range_shell[1]*SHELL_FACTOR)

    # objective
    model.add_objective('wing.weight', scaler=1)

    # constraint
    model.add_constraint('wing.stress', lower=0., upper=max_shear_strength * STRESS_FAC)
    model.add_subsystem('con_cmp1', ExecComp('con1 = (ribs * '+str(RIB_FACTOR)+') - int(ribs[0] * '+str(RIB_FACTOR)+')'))
    model.add_constraint('con_cmp1.con1', upper=0.1)

    prob = Problem(model)

    # setup the optimization
    #prob.driver = pyOptSparseDriver() #ScipyOptimizeDriver()
    #prob.driver.options['optimizer'] = 'SLSQP' #'SLSQP'
    #prob.driver.options['tol'] = 1e-6
    #prob.driver.opt_settings = {'eps': 1e-6}
    #prob.driver.options['maxiter'] = 100000
    #prob.driver.options['disp'] = True


    prob.driver =  ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'  # ['Nelder-Mead', 'Powell', 'CG', 'BFGS', 'Newton-CG', 'L-BFGS-B', 'TNC', 'COBYLA', 'SLSQP']
    prob.driver.options['tol'] = 1e-6
    prob.driver.opt_settings = {'eps': 1e-6}
    prob.driver.options['maxiter'] = 1000
    prob.driver.options['disp'] = True


    prob.setup()
    prob.set_solver_print(level=0)
    #prob.model.approx_totals()
    prob.setup(check=True, mode='fwd')
    prob.run_driver()



    '''
    # setup recorder
    recorder = SqliteRecorder(Constants().WORKING_DIR + '/openMdaoLog.sql')
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options['record_desvars'] = True
    prob.driver.recording_options['record_responses'] = True
    prob.driver.recording_options['record_objectives'] = True
    prob.driver.recording_options['record_constraints'] = True
    '''

    print('done')
    print('ribs: ' + str(prob['wing.ribs'] / RIB_FACTOR))
    print('shell: ' + str(prob['wing.shell'] / SHELL_FACTOR) + ' m')

    print('weight= ' + str(prob['wing.weight'] / WEIGHT_FAC))
    print('stress= ' + str(prob['wing.stress'] / STRESS_FAC))

    print('execution counts wing: ' + str(prob.model.wing.executionCounter))


if __name__ == '__main__':
    run_open_mdao()