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
from scipy import optimize

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/OpenMDAO')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyDOE2')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pyoptsparse')
from openmdao.api import Problem, ExecComp, pyOptSparseDriver, ScipyOptimizeDriver, IndepVarComp, ExplicitComponent, SqliteRecorder, ScipyKrylov, Group, DirectSolver, NewtonSolver, NonlinearBlockGS
from openmdao.core.problem import Problem
from openmdao.core.indepvarcomp import IndepVarComp

from wingConstruction.wingUtils.constants import Constants
from wingConstruction.multi_run import MultiRun
from wingConstruction.wingUtils.defines import *
from myutils.plot_helper import PlotHelper
from myutils.time_track import TimeTrack

PROJECT_NAME_PREFIX = 'newtonOpti'

LOG_FILE_PATH = Constants().WORKING_DIR + '/' + PROJECT_NAME_PREFIX + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
SHELL_FACTOR = 1e-2
RIB_FACTOR = 1e-6
WEIGHT_FAC = 1e-3
STRESS_FAC = 1e-8

WEIGHT_PANALTY_FAC = 0

USE_ABA = True

PGF = True


class NewtonOpt:

    def __init__(self):
        pass

    def opti_it(self, rib_range=[16,17,18,19,20,21,22,23]):
        ######################
        ### needed Objects ###
        self.runner = MultiRun(use_calcu=not USE_ABA, use_aba=USE_ABA, non_liner=False, project_name_prefix=PROJECT_NAME_PREFIX, force_recalc=False)
        self.executionCounter = 0
        self.timer = TimeTrack()

        opti_ribs = []
        opti_shell = []
        opti_stress = []
        opti_weights = []
        for r in rib_range:
            r = int(r)
            init_guess = (range_shell[0] + range_shell[1]) / 2.
            self.executionCounter = 0
            root = optimize.newton(self.shell_predict, init_guess, args=[r], tol=1.48e-08, maxiter=50)
            opti_ribs.append(r)
            opti_shell.append(root)
            stress, weight = self.calc_stress_weight(r, root)
            opti_stress.append(stress)
            opti_weights.append(weight)
            print('execution count: {:d}'.format(self.executionCounter))
        print('DONE')
        print(opti_ribs)
        print(opti_shell)
        print(opti_stress)
        print(opti_weights)
        best_i = opti_weights.index(min(opti_weights))
        print('BEST:')
        print('ribs:' + str(opti_ribs[best_i]))
        print('shell:' + str(opti_shell[best_i]))
        print('stress:' + str(opti_stress[best_i]))
        print('weight:' + str(opti_weights[best_i]))

    def shell_predict(self, shell_thick, rib_num):
        stress, _ = self.calc_stress_weight(shell_thick, rib_num)
        self.executionCounter += 1
        return stress - max_shear_strength

    def calc_stress_weight(self, shell_thick, rib_num):
        self.runner.project_name_prefix = PROJECT_NAME_PREFIX + '_{:05d}'.format(self.executionCounter)
        pro = self.runner.new_project_r_t(rib_num, shell_thick)
        pro = self.runner.run_project(pro, used_cpus=1)
        res = pro.resultsCalcu
        if USE_ABA:
            res = pro.resultsAba
        stress = res.stressMisesMax * STRESS_FAC
        weight = pro.calc_wight() * WEIGHT_FAC
        return stress, weight


if __name__ == '__main__':
    nw = NewtonOpt()
    nw.opti_it()