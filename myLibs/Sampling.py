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

import matplotlib.pyplot as plt
import numpy as np
import math
import sys


class Sampling:

    def __init__(self):
        print('done')

    def bool_mat_to_list(self, mat):
        samples = []
        for r in range(0, mat.shape[0]):
            for c in range(0, mat.shape[0]):
                if mat[c][r] == 1:
                    samples.append([c, r])
        return samples

    '''
    blocks x dim design (with n x m sample variable steps)
    it has to be a n x n design space
    '''
    def enhanced_latin_hypercube(self, n):
        cube_size = n
        if math.sqrt(cube_size) > 0.:
            cube_size = (int(math.sqrt(cube_size))+1)**2
        samples = self.enhanced_latin_hypercube_n_x_n(cube_size)
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
            print('kicked out: {:d}, {:d}'.format(max_dist_pos[0], max_dist_pos[1]))
        return samples

    '''
    generates enhanced latin hypercube design space
    represented as n x n matrix (filled with 0)
    where the sample places are marked with 1
    it has to be a n x n design space where n is int(x)**2
    '''
    def enhanced_latin_hypercube_n_x_n(self, n):
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


if __name__ == '__main__':
    sam = Sampling()
    sample_mat = sam.enhanced_latin_hypercube(14)
    xy = sam.bool_mat_to_list(sample_mat)
    plotXY = np.array(xy).T.tolist()
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(plotXY[0], plotXY[1], 'o-')
    ax.set_xticks(range(0,14), minor=False)
    ax.set_yticks(range(0,14), minor=False)
    ax.grid(True)
    plt.show()