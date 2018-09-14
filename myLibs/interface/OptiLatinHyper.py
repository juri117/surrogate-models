__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Sampling plans
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import numpy as np
from ast import literal_eval
import math
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../../lib/pyKriging')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../../lib/inspyred')

import pyKriging
from pyKriging.krige import kriging
from pyKriging.samplingplan import samplingplan

from utils.PlotHelper import PlotHelper


class OptiLatinHyper:

    def __init__(self):
        pass

    def generate_sample_plan(self, point_count, dimension, bounds, base=None):
        sp = samplingplan(dimension)
        X = sp.optimallhc(point_count)
        norm_point = X
        points = []
        for i in range(0, point_count):
            scaled_point = []
            for d in range(0, dimension):
                scaled_point.append(bounds[d][0] + (norm_point[i][d] * (bounds[d][1] - bounds[d][0])))
            points.append(scaled_point)
        return points


if __name__ == '__main__':
    sam = OptiLatinHyper()

    import matplotlib.pyplot as plt

    samples = sam.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)])
    for i in range(0, 14):
        plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    plt.show()