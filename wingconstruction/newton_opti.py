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

from wingconstruction.wingutils.constants import Constants
from wingconstruction.multi_run import MultiRun
from wingconstruction.wingutils.defines import *
from myutils.plot_helper import PlotHelper
from myutils.time_track import TimeTrack

PROJECT_NAME_PREFIX = 'newtonOpti'

LOG_FILE_PATH = Constants().WORKING_DIR + '/' + PROJECT_NAME_PREFIX + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'


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
        self.write_newton_log('iter,time,ribs,shell,stress,weight')
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
            stress, weight = self.calc_stress_weight(root, r)
            opti_stress.append(stress)
            opti_weights.append(weight)
            print('execution count: {:d}'.format(self.executionCounter))
            self.write_newton_log(str(self.executionCounter) + ','
                           + str(self.timer.get_time()) + ','
                           + str(r) + ','
                           + str(root) + ','
                           + str(stress) + ','
                           + str(weight))
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
        if shell_thick > range_shell[1]:
            shell_thick = range_shell[1]
        if shell_thick < range_shell[0]:
            shell_thick = range_shell[0]
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
        stress = res.stressMisesMax
        weight = pro.calc_wight()
        return stress, weight

    def write_newton_log(out_str):
        out_str = out_str.replace('[', '')
        out_str = out_str.replace(']', '')
        output_f = open(LOG_FILE_PATH, 'a')  # 'a' so we append the file
        output_f.write(out_str + '\n')
        output_f.close()


if __name__ == '__main__':
    nw = NewtonOpt()
    nw.opti_it()