__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional DoE
# date            :2018-079-20
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import numpy as np
from ast import literal_eval
import math
import sys
from utils.PlotHelper import PlotHelper


class DoE:

    def __init__(self, input_names, ranges, func, level_count=2):
        self._inputNames = input_names
        self._ranges = ranges
        self._func = func
        self._k = len(input_names)
        self._l = level_count
        self._levels = self._ranges
        self._results = []
        self.auto_define_levels()

    def auto_define_levels(self):
        levels = []
        for i in range(0, self._k):
            inc = (self._ranges[i][1] - self._ranges[i][0]) / (self._l - 1)
            lev = []
            for l in range(0, self._l):
                lev.append(self._ranges[i][0] + (l * inc))
            levels.append(lev)
        self._levels = levels

    def corellation(self):
        inp = np.zeros(self._k)
        for i in range(0, self._l**self._k):
            bin = self.base(i, self._l)
            while len(bin) < self._k:
                bin = [0] + bin
            self._run(bin)

    def print_res_table(self):
        ref = self._results[0].res
        out = ''
        for i in range(self._k):
            out += '{:s}\t\t|'.format(self._inputNames[i])
        out += '{:s}\t\t\t\t|{:s}\n'.format('result', '%')
        for r in self._results:
            for i in range(self._k):
                out += '{:f}\t|'.format(self._levels[i][r.inp[i]])
            out += '{:f}\t|{:f}\n'.format(r.res, 100*(r.res/ref))
        print(out)

    '''
    :param decimal int number
    :param base the base to use for convertion
    :return the decimal number noted as base-number
    '''
    def base(self, decimal, base):
        other_base = []
        while decimal != 0:
            other_base.append(int(decimal % base))
            decimal = decimal / base
        out = []
        for s in reversed(other_base):
            if s != 0 or len(out) > 0:
                out.append(s)
        return out

    def _run(self, inp):
        inputs = []
        out = ''
        for i in range(0, self._k):
            if inp[i] == 1:
                out += self._inputNames[i]
            inputs.append(self._levels[i][inp[i]])
        print('['+out+']')
        fx = self._func(inputs)
        res = DoERun(inp, fx)
        self._results.append(res)


class DoERun:
    def __init__(self, inp, res):
        self.inp = inp
        self.res = res


if __name__ == '__main__':
    input_names = ['ribs', 'shell']
    from wingConstruction.utils.defines import *
    ranges = [range_rib, range_shell]
    from wingConstruction.MultiRun import MultiRun
    runner = MultiRun(use_calcu=True, use_aba=False, non_liner=False, force_recalc=False, project_name_prefix='DoE')
    d = DoE(input_names, ranges, runner.calc_stress, level_count=3)
    d.corellation()
    d.print_res_table()
