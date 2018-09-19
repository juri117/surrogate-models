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

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/OpenMDAO')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyDOE2')
from openmdao.api import Problem, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingConstruction.utils.Constants import Constants
from wingConstruction.MultiRun import MultiRun
from wingConstruction.utils.defines import *

LOG_FILE_PATH = Constants().WORKING_DIR + '/om_iterations_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'

PROJECT_NAME_PREFIX = 'iter'

SHELL_FACTOR = 1
RIB_FACTOR = 1e-6

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
        self.add_output('stress', val=.2)
        self.add_output('weight', val=.2)

        self.declare_partials('*', '*', method='fd')
        self.executionCounter = 0

    def compute(self, inputs, outputs):
        ribs = int(round(inputs['ribs'][0] / RIB_FACTOR))
        shell = inputs['shell'][0] / SHELL_FACTOR
        pro = self.runner.new_project_r_t(ribs, shell)
        pro = self.runner.run_project(pro)
        outputs['stress'] = pro.resultsCalcu.stressMisesMax
        outputs['weight'] = pro.calc_wight()
        write_to_log(str(self.executionCounter) + ','
                     + datetime.now().strftime('%H:%M:%S') + ','
                     + str(inputs['ribs']) + ','
                     + str(inputs['shell']) + ','
                     + str(outputs['stress']) + ','
                     + str(outputs['weight']))
        self.executionCounter += 1
        print('{:f}, {:f} -> {:f}, {:f}'.format(inputs['ribs'][0], inputs['shell'][0], outputs['stress'][0], outputs['weight'][0]))



def write_to_log(outStr):
    outStr = outStr.replace('[', '')
    outStr = outStr.replace(']', '')
    outputF = open(LOG_FILE_PATH, 'a') #'a' so we append the file
    outputF.write(outStr + '\n')
    outputF.close()

def runOpenMdao():
    prob = Problem()

    indeps = prob.model.add_subsystem('indeps', IndepVarComp(), promotes=['*'])
    indeps.add_output('ribs', range_rib[1] * RIB_FACTOR)
    indeps.add_output('shell', (range_shell[1])*SHELL_FACTOR)

    prob.model.add_subsystem('wing', WingStructure())
    prob.model.connect('ribs', 'wing.ribs')
    prob.model.connect('shell', 'wing.shell')

    # setup the optimization
    prob.driver = ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['tol'] = 1e-6
    prob.driver.options['maxiter'] = 100000

    # setup recorder
    recorder = SqliteRecorder(Constants().WORKING_DIR + '/openMdaoLog.sql')
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options['record_desvars'] = True
    prob.driver.recording_options['record_responses'] = True
    prob.driver.recording_options['record_objectives'] = True
    prob.driver.recording_options['record_constraints'] = True

    # design variables, limits and constraints
    prob.model.add_design_var('ribs', lower=range_rib[0]*RIB_FACTOR, upper=range_rib[1]*RIB_FACTOR)
    prob.model.add_design_var('shell', lower=range_shell[0]*SHELL_FACTOR, upper=range_shell[1]*SHELL_FACTOR)

    # objective
    prob.model.add_objective('wing.weight', scaler=1)

    # constraint
    prob.model.add_constraint('wing.stress', lower=0., upper=max_shear_strength)

    prob.setup()
    prob.set_solver_print(level=0)
    prob.model.approx_totals()
    prob.run_driver()

    print('done')
    print('ribs: ' + str(prob['wing.ribs']))
    print('cabin angle: ' + str(prob['wing.shell']) + ' m')

    print('weight= ' + str(prob['wing.weight']))
    print('stress= ' + str(prob['wing.stress']))

    print('execution counts wing: ' + str(prob.model.wing.executionCounter))

if __name__ == '__main__':
    runOpenMdao()