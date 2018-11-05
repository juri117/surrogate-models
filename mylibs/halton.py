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


class Halton:

    def __init__(self):
        pass

    def prime(self, n):
        """
        :param n
        :return the n-th prime number (starting from n=0 => 2)
        """
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

    def halton(self, i, prim):
        """
        :param i the index of the generated number
        :param prim the base
        """
        # add one to exclude 0 as result
        i = int(i+1)
        baseStr = self.base(i, prim)
        baseStrInv = self.str_inverse(baseStr)
        decimal = self.base_fract_str_to_float(baseStrInv, prim)
        return decimal

    def str_inverse(self, nums):
        """
        :param number list
        :return the inverse list of numbers
        """
        str_inv = np.zeros((len(nums)))
        for i in range(0, len(nums)):
            str_inv[len(nums)-1 - i ] = nums[i]
        return list(str_inv)

    def base_fract_str_to_float(self, base_fracs, base):
        """
        :param baseFracs a list that represents a number with base notation
        :param the base of baseFracStr
        :return a float of baseFracs evaluated with base, as fraction
        """
        res = 0.
        for i in range(0, len(base_fracs)):
            res += int(base_fracs[i]) * int(base) ** (-1 * (i + 1))
        return res

    def base(self, decimal, base):
        """
        :param decimal int number
        :param base the base to use for convertion
        :return the decimal number noted as base-number
        """
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
        """
        generates sampling plan
        :param point_count: number of sampling points
        :param dimension: dimension of the sampling plan
        :param bounds: vector of tooples representing the bounds for every input
        :param base: vector of bases for every input
        :return: matrix: list of point_count entries with each dimension entries representing the sampling plan
        """
        if not base is None:
            used_base = base
        else:
            used_base = []
            for i in range(0, dimension):
                used_base.append(self.prime(i))
        points = []
        for i in range(0, point_count):
            scaled_point = []
            for d in range(0, dimension):
                halt = self.halton(i, used_base[d])
                scaled_point.append(bounds[d][0] + (halt * (bounds[d][1]- bounds[d][0])))
            points.append(scaled_point)
        return points
