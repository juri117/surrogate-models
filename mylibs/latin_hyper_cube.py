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
from myutils.plot_helper import PlotHelper
import operator


class LatinHyperCube:

    def __init__(self):
        pass

    '''
    def bool_mat_to_list(self, mat):
        samples = []
        for r in range(0, mat.shape[0]):
            for c in range(0, mat.shape[0]):
                if mat[c][r] == 1:
                    samples.append([r, c])
        return samples
    '''

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


    '''
    def enhanced_latin_hypercube(self, n):
        cube_size = n
        if math.sqrt(cube_size) % 1 > 0.:
            cube_size = (int(math.sqrt(cube_size))+1)**2
        samples = self.enhanced_latin_hypercube_k_pow_x(2, cube_size)
        #i_far_away = [0, 0]
        #dist_at_i = -1.
        center = (cube_size-1)/2
        dist_mat = np.zeros((cube_size, cube_size))
        for r in range(0, samples.shape[0]):
            for c in range(0, samples.shape[1]):
                if samples[c][r] == 1:
                    dist_mat[c][r] = math.sqrt((r-center)**2 + (c-center)**2)
        while samples.shape[0] > n:
            max_dist_pos = np.unravel_index(dist_mat.argmax(), dist_mat.shape)
            samples = np.delete(samples, max_dist_pos[0], axis=0)
            samples = np.delete(samples, max_dist_pos[1], axis=1)
            dist_mat = np.delete(dist_mat, max_dist_pos[0], axis=0)
            dist_mat = np.delete(dist_mat, max_dist_pos[1], axis=1)
            #print('kicked out: {:d}, {:d}'.format(max_dist_pos[1], max_dist_pos[0]))
        return samples
    '''

    def calc_dist(self, coords, target):
        sum = 0
        for c in coords:
            sum += (c - target)**2
        return math.sqrt(sum)

    '''
    blocks x dim design (with n x m sample variable steps)
    it has to be a n x n design space
    '''
    def enhanced_latin_hypercube(self, k, n):
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


    '''
    def enhanced_latin_hypercube_2_pow_x(self, n):
        edge_devision = math.sqrt(n)
        if edge_devision % 1 > 0.:
            print('ERROR, n has to be x^2')
            sys.exit(0)
        edge_devision = int(edge_devision)
        samples = np.zeros((n, n))
        row_count = 0
        col_count = 0
        for i in range(0, edge_devision):
            for e in range(0, edge_devision):
                #samples.append([(e*4) + col_count, row_count])
                samples[row_count][(e*edge_devision) + col_count] = 1
                row_count += 1
            col_count += 1
        return samples
    '''

    '''
    generates enhanced latin hypercube design space
    represented as n x n matrix (filled with 0)
    where the sample places are marked with 1
    it has to be a n x n design space where n is int(x)**2
    '''
    def enhanced_latin_hypercube_k_pow_x(self, k, n):
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

    def generate_sample_plan(self, point_count, dimension, bounds, base=None):
        #if dimension != 2:
        #    print('LatinHyperCube does not support dimensions other than 2 yet!')
        #    sys.exit(9061548)
        sample_mat = self.enhanced_latin_hypercube(dimension, point_count)
        #sample_mat = self.enhanced_latin_hypercube_rn(dimension, point_count)
        sample_indices = self.bool_mat_to_list(sample_mat)
        normPoint = list(np.array(sample_indices) * (1/(point_count-1)))
        points = []
        for i in range(0, point_count):
            scaledPoint = []
            for d in range(0, dimension):
                scaledPoint.append(bounds[d][0] + (normPoint[i][d] * (bounds[d][1] - bounds[d][0])))
            points.append(scaledPoint)
        return points

if __name__ == '__main__':
    #import matplotlib.pyplot as plt
    sam = LatinHyperCube()

    import matplotlib.pyplot as plt
    samples = sam.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05), (1, 5)])
    #samples = sam.enhanced_latin_hypercube_2_pow_x(16)
    for i in range(0, len(samples)):
        plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    plt.show()

    pltHalton = PlotHelper([], fancy=True, pgf=False)
    import matplotlib.pyplot as plt

    ax1 = pltHalton.fig.add_subplot(121)
    ax2 = pltHalton.fig.add_subplot(122)

    sample_mat_full = sam.enhanced_latin_hypercube(2, 16)
    xy_full = sam.bool_mat_to_list(sample_mat_full)
    plotXY_full = np.array(xy_full).T.tolist()
    ax1.plot(plotXY_full[0], plotXY_full[1], 'bo', markersize=5)

    deleteX = [0, 15, 4, 1]
    deleteY = [0, 15, 1, 4]
    for i in range(0, len(deleteX)):
        #ax1.plot([deleteX[i]], [deleteY[i]], 'rx', mew=2, ms=10)
        ax1.plot([-1, 16], [deleteY[i], deleteY[i]], 'r-', linewidth=2)
        ax1.plot([deleteX[i], deleteX[i]], [-1, 16], 'r-', linewidth=2)
        ax1.text(deleteX[i]-0.3, deleteY[i]-0.3, str(i+1), pltHalton.font, fontweight='bold')

    ax1.set_xticks(range(0, 16), minor=False)
    ax1.set_yticks(range(0, 16), minor=False)
    ax1.xaxis.set_ticklabels([])
    ax1.yaxis.set_ticklabels([])
    ax1.set_xlim(16, -1)
    ax1.set_ylim(16, -1)
    ax1.grid(True)

    #ax = fig.add_subplot(1, 2, 2)

    sample_mat = sam.enhanced_latin_hypercube(2, 12)
    xy = sam.bool_mat_to_list(sample_mat)
    plotXY = np.array(xy).T.tolist()

    ax2.plot(plotXY[0], plotXY[1], 'bo', markersize=5)
    ax2.set_xticks(range(0,12), minor=False)
    ax2.set_yticks(range(0,12), minor=False)
    ax2.xaxis.set_ticklabels([])
    ax2.yaxis.set_ticklabels([])
    ax2.set_xlim(12, -1)
    ax2.set_ylim(12, -1)
    ax2.grid(True)
    ax2.grid(True)
    #plt.show()
    pltHalton.fig.set_size_inches(5, 2.5)
    plt.tight_layout()
    pltHalton.save('../data_out/plot/latinHyper.pdf')
    pltHalton.show()