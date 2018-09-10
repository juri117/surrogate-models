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
from utils.PlotHelper import PlotHelper


class Halton:

    def __init__(self):
        print('done')

    '''
    :param n
    :return the n-th prime number (starting from n=0 => 2)
    '''
    def prime(self, n):
        n += 1
        prime_list = [2]
        num = 3
        while len(prime_list) < n:
            for p in prime_list:
                if num % p == 0:
                    break
            else:
                prime_list.append(num)
            num += 2
        return prime_list[-1]

    '''
    :param i the index of the generated number
    :param prim the base
    '''
    def halton(self, i, prim):
        # add one to exclude 0 as result
        i = int(i+1)
        baseStr = self.base(i, prim)
        baseStrInv = self.str_inverse(baseStr)
        decimal = self.base_fract_str_to_float(baseStrInv, prim)
        print(decimal)
        return decimal

    '''
    :param number list
    :return the inverse list of numbers
    '''
    def str_inverse(self, nums):
        strInv = np.zeros((len(nums)))
        for i in range(0, len(nums)):
            strInv[len(nums)-1 - i ] = nums[i]
        return list(strInv)

    '''
    :param baseFracs a list that represents a number with base notation
    :param the base of baseFracStr
    :return a float of baseFracs evaluated with base, as fraction
    '''
    def base_fract_str_to_float(self, baseFracs, base):
        res = 0.
        for i in range(0, len(baseFracs)):
            res += int(baseFracs[i]) * int(base) ** (-1 * (i + 1))
        return res

    '''
    :param decimal int number
    :param base the base to use for convertion
    :return the decimal number noted as base-number
    '''
    def base(self, decimal, base):
        #list = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        other_base = []
        while decimal != 0:
            other_base.append(int(decimal % base))
            decimal = decimal / base
        out = []
        for s in reversed(other_base):
            if s != 0 or len(out) > 0:
                out.append(s)
        return out


    def generate_sample_plan(self, point_count, dimension, bounds, base=None):
        if not base is None:
            usedBase = base
        else:
            usedBase = []
            for i in range(0, dimension):
                usedBase.append(self.prime(i))
        points = []
        for i in range(0, point_count):
            scaledPoint = []
            for d in range(0, dimension):
                halt = self.halton(i, usedBase[d])
                scaledPoint.append(bounds[d][0] + (halt * (bounds[d][1]- bounds[d][0])))
            points.append(scaledPoint)
        return points



if __name__ == '__main__':
    hal = Halton()
    print(hal.halton(10, 2))

    #for i in range(0, 10):
    #    point = hal.halton(i, 2)

    #sys.exit(15)

    pltHalton = PlotHelper([], fancy=True, pgf=True)
    import matplotlib.pyplot as plt

    ax1 = pltHalton.fig.add_subplot(121)
    ax2 = pltHalton.fig.add_subplot(122)

    #ax = fig.add_subplot(1, 2, 1)
    for i in range(0, 100):
        point = [hal.halton(i, 11), hal.halton(i, 29)]
        ax1.plot([point[1]], [point[0]], 'bo', markersize=3)
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])

    #ax = fig.add_subplot(1, 2, 2)
    for i in range(0, 100):
        point = [hal.halton(i, 2), hal.halton(i, 19)]
        ax2.plot([point[1]], [point[0]], 'bo', markersize=3)
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticklabels([])

    pltHalton.fig.set_size_inches(5, 2.5)
    plt.tight_layout()
    pltHalton.save('../dataOut/halton.pdf')
    pltHalton.show()

    #samples = hal.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)])
    #for i in range(0, 14):
    #    plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    #plt.show()