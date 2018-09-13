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
                scaled_point.append(bounds[d][0] + (norm_point[i][d] * (bounds[d][1] - bounds[d][0])))
            points.append(scaled_point)
        return points


if __name__ == '__main__':
    sam = StructuredSample()

    import matplotlib.pyplot as plt

    samples = sam.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)])
    for i in range(0, 14):
        plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    plt.show()
