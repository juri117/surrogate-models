__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :validation of surrogates
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

import matplotlib.pyplot as plt
import numpy as np
import math
import sys


class Validation:

    def __init__(self):
        print('done')

    '''
    calc the deviation of a matrix with known solutions to the surrogate solution
    '''
    def calc_deviation(self, x, y, values, surro_func):
        count = 0
        sum_deviation = 0
        # sample_indices = np.array([known_x_i, known_y_i]).T.tolist()
        for i_r in range(0, len(x)):
            for i_s in range(0, len(y)):
                devi = values[i_s][i_r] - surro_func([x[i_r], y[i_s]])
                sum_deviation += abs(devi)
                count += 1
        avg_deviation = sum_deviation / count
        avg_deviation_per = avg_deviation / np.array(values).mean()
        print('avg deviation: {:.3e} (-> {:.3f}%)'.format(avg_deviation, avg_deviation_per * 100.))
        return avg_deviation_per


if __name__ == '__main__':
    val = Validation()
