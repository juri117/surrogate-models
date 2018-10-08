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
from myUtils.PlotHelper import PlotHelper


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

    str_plt = PlotHelper(['', ''], fancy=True, pgf=True)
    import matplotlib.pyplot as plt

    samples = sam.generate_sample_plan(16, 2, [(0, 30), (0, 30)])
    for i in range(0, 16):
        str_plt.ax.plot([samples[i][0]], [samples[i][1]], 'bo')

    str_plt.ax.xaxis.set_ticklabels([])
    str_plt.ax.yaxis.set_ticklabels([])
    #str_plt.ax.set_xticks(range(0, 31), minor=False)
    #str_plt.ax.set_yticks(range(0, 31), minor=False)
    str_plt.ax.locator_params(nbins=4, axis='x')
    str_plt.ax.locator_params(nbins=4, axis='y')
    str_plt.ax.grid(True)

    str_plt.finalize(width=1.5, height=1.5, show_legend=False)
    str_plt.save('../dataOut/plot/structuredSamp.pdf')
    str_plt.show()
