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
import math
import sys


class LatinHyperCube:

    def __init__(self):
        pass

    def bool_mat_to_list(self, mat):
        n = mat.shape[0]
        k = len(mat.shape)
        samples = []
        i_s = np.zeros(k)
        while True:
            i_s = [int(x) for x in i_s]
            if mat[tuple(i_s)] == 1:
                samples.append(np.flip(np.array(i_s), 0))
            if not any(x < n - 1 for x in i_s):
                break
            self.increase_edge(i_s, n)
        return samples

    def calc_dist(self, coords, target):
        sum = 0
        for c in coords:
            sum += (c - target)**2
        return math.sqrt(sum)

    def enhanced_latin_hypercube(self, k, n):
        """
        blocks x dim design (with n x m sample variable steps)
        it has to be a n x n design space
        :param k: dimension
        :param n: number of points
        :return: list of sampling points
        """
        cube_size = n
        if n**(1./k) % 1 > 0.:
            cube_size = (int(n**(1./k)) + 1) ** k
        samples = self.enhanced_latin_hypercube_k_pow_x(k, cube_size)
        dist_mat = samples.copy()
        center = (cube_size - 1) / 2

        i_s = np.zeros(k)
        while True:
            i_s = [int(x) for x in i_s]
            if dist_mat[tuple(np.flip(i_s, 0))] == 1:
                dist_mat[tuple(np.flip(i_s, 0))] = self.calc_dist(np.flip(i_s, 0), center)
            if not any(x < cube_size - 1 for x in i_s):
                break
            self.increase_edge(i_s, cube_size)
        while samples.shape[0] > n:
            max_dist_pos = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
            for d in range(0, k):
                samples = np.delete(samples, max_dist_pos[d], axis=d)
                dist_mat = np.delete(dist_mat, max_dist_pos[d], axis=d)
        return samples


    def enhanced_latin_hypercube_k_pow_x(self, k, n):
        """
        generates enhanced latin hypercube design space
        represented as n x n matrix (filled with 0)
        where the sample places are marked with 1
        it has to be a n x n design space where n is int(x)**2
        :param k: dimension
        :param n: number of points
        :return: list of sampling points
        """
        edge_devision = n**(1./k)
        if edge_devision % 1 > 0.000000001 and edge_devision % 1 < 0.999999999:
            print('ERROR, n has to be x^k')
            sys.exit(0)
        edge_devision = int(round(edge_devision))
        edge_devision = int(edge_devision)
        dimensions = ()
        for ik in range(0, k):
            dimensions += (n,)
        samples = np.zeros(dimensions)
        cube_is = np.zeros(k)
        edge_is = np.zeros(k)
        while True:
            i_s = np.flip(edge_is, 0) + (cube_is * edge_devision)
            i_s = [int(x) for x in i_s]
            samples[tuple(i_s)] = 1
            if not any(x < edge_devision-1 for x in edge_is):
                break
            self.increase_edge(cube_is, edge_devision)
            self.increase_edge(edge_is, edge_devision)
        return samples

    def increase_edge(self, edge_count, edge):
        for i in range(0, len(edge_count)):
            # -1 here because we start indexing at 0
            if edge_count[i] >= edge - 1:
                edge_count[i] = 0
            else:
                edge_count[i] += 1
                return edge_count
        return edge_count

    def generate_sample_plan(self, point_count, dimension, bounds):
        """
        generates sampling plan
        :param point_count: number of sampling points
        :param dimension: dimension of the sampling plan
        :param bounds: vector of tooples representing the bounds for every input
        :return: matrix: list of point_count entries with each dimension entries representing the sampling plan
        """
        sample_mat = self.enhanced_latin_hypercube(dimension, point_count)
        sample_indices = self.bool_mat_to_list(sample_mat)
        norm_point = list(np.array(sample_indices) * (1/(point_count-1)))
        points = []
        for i in range(0, point_count):
            scaled_point = []
            for d in range(0, dimension):
                scaled_point.append(bounds[d][0] + (norm_point[i][d] * (bounds[d][1] - bounds[d][0])))
            points.append(scaled_point)
        return points
