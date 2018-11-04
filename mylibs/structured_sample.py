__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Sampling plan full-factorial sampling
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import numpy as np
import math


class StructuredSample:

    def __init__(self):
        pass

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
        sample_indices = []
        edge = math.ceil(point_count**(1/float(dimension)))
        edge_count = np.zeros((dimension))

        while len(sample_indices) < point_count:
            p = []
            for di in range(0, dimension):
                p.append(edge_count[di])
            sample_indices.append(p)
            edge_count = self.increase_edge(edge_count, edge)

        norm_point = list(np.array(sample_indices) * (1 / (edge-1)))
        points = []
        for i in range(0, point_count):
            scaled_point = []
            for d in range(0, dimension):
                #fix for if one row is missing
                point = norm_point[i][d]
                if max(np.array(norm_point)[:,d]) < 1 and point_count > 2:
                    point = point * (1. / max(np.array(norm_point)[:,d]))

                scaled_point.append(bounds[d][0] + (point * (bounds[d][1] - bounds[d][0])))
            points.append(scaled_point)
        return points
