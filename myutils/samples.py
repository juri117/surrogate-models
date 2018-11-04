__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :definition of sample functions for testing
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import math


def f_3d(coords):
    x = coords[0]
    y = coords[1]
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4 + 0.05*y**2 - 0.001*y**4 - 0.005*y**3


def f_2d(x):
    return 6 - (math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4)


def f_2d_big(x):
    return math.sin(x)+3.+0.035*x**2+0.003*x**3-0.00045*x**4+0.0000005*x**6


def generate_sample_data(f_3d, px_edge, py_edge):
    px = []
    py = []
    pz = []
    for iX in range(0, len(px_edge)):
        for iY in range(0, len(py_edge)):
            px.append(px_edge[iX])
            py.append(py_edge[iY])
            pz.append(f_3d([px_edge[iX], py_edge[iY]]))
    return px, py, pz


def print_mat(mat, octave=False):
    if octave:
        print('['+';\n'.join([''.join(['{:10.3} '.format(item) for item in row]) for row in mat])+']')
    else:
        print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))
