
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import scipy

from RBF import RBF
from utils.PlotHelper import PlotHelper
from utils.samples import *

FONT_SIZE=14

def printMat(mat):
    print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4

# gaus RBF (GS)
def gausRBF(a, r):
    return math.e**(-((a*r)**2))

def multiQuadRBF(a, r):
    return math.sqrt(1 + (a * r) ** 2)

def invMultiQuadRBF(a, r):
    return (1+r**2)**(a/2)

# todo: other RBFunctions here...

def rbf_calc_coefficiants(rbf_const, knownX, knownVal, rbf):
    knownCount = len(knownX)
    mat = np.zeros((knownCount, knownCount))
    radialMat = np.zeros((knownCount, knownCount))
    for i1 in range(0, knownCount):
        for i2 in range(i1, knownCount):
            radius = math.sqrt((knownX[i1] - knownX[i2]) ** 2.)
            radialMat[i1][i2] = radius
            mat[i1][i2] = rbf(rbf_const, radius)
            mat[i2][i1] = mat[i1][i2]
    print('radialMat:')
    printMat(radialMat)
    print('mat:')
    printMat(mat)
    return np.linalg.solve(mat, knownVal)

def rbfSolution(x, knownX, coefficients, rbf, rbf_const):
    n = len(knownX)
    res = 0.
    for i in range(0, n):
        radius = math.sqrt((x - knownX[i])**2)
        res += coefficients[i] * rbf(rbf_const, radius)
    return res

if __name__ == '__main__':
    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=False)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f,fx))

    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # now we pretend we only know a view points
    knownParams = [0., 2., 4., 6., 8., 10.]
    knownValues = list(map(f, knownParams))

    fx = fx.reshape((len(fx), 1))
    # RBF 1
    rbf1 = RBF([knownParams], knownValues)
    a1 = 1.
    rbf1.update_param(a1, 'gaus')
    rbfY1 = list(map(rbf1.predict, fx))
    plt1.ax.plot(fx, rbfY1, 'b--', label=r'$f_{gaus-RBF}$ mit $a = ' + str(a1) + '$')

    # RBF 2
    rbf2 = RBF([knownParams], knownValues)
    a2 = 0.24
    rbf2.update_param(a2, 'gaus')
    rbfY2 = list(map(rbf2.predict, fx))
    plt1.ax.plot(fx, rbfY2, 'b:', label=r'$f_{gaus-RBF}$ mit $a = ' + str(a2) + '$')

    # RBF 3
    rbf3 = RBF([knownParams], knownValues)
    a3 = 1.
    rbf3.update_param(a3, 'inverse-multi-quadratic')
    rbfY3 = list(map(rbf3.predict, fx))
    plt1.ax.plot(fx, rbfY3, 'c--', label=r'$f_{imq-RBF}$ mit $a = ' + str(a3) + '$')

    # RBF 4
    rbf4 = RBF([knownParams], knownValues)
    a4 = 1.8
    rbf4.update_param(a4, 'inverse-multi-quadratic')
    rbfY4 = list(map(rbf4.predict, fx))
    plt1.ax.plot(fx, rbfY4, 'c:', label=r'$f_{imq-RBF}$ mit $a = ' + str(a4) + '$')

    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    plt1.finalize(width=8, height=5, legendLoc='lower left', legendNcol=2)
    plt1.ax.legend(loc='lower left', ncol=2)  # , mode="expand")
    #plt1.save('dataOut/radialBasisR2.svg')
    #plt1.save('dataOut/radialBasisR2.pdf')
    plt1.show()


    #a1 = 1.8
    #rbf1 = invMultiQuadRBF
    #a1 = 1.
    #rbf1 = gausRBF
    #coeff1 = rbf_calc_coefficiants(a1, knownParams, knownValues, rbf1)

    #a2 = 0.24
    #rbf2 = gausRBF
    #coeff2 = rbf_calc_coefficiants(a2, knownParams, knownValues, rbf2)

    #a3 = 1.
    #rbf3 = invMultiQuadRBF
    #coeff3 = rbf_calc_coefficiants(a3, knownParams, knownValues, rbf3)

    #a4 = 1.8
    #rbf4 = invMultiQuadRBF
    #coeff4 = rbf_calc_coefficiants(a4, knownParams, knownValues, rbf4)

    #rbfy1 = np.zeros((len(fx), 1))
    #rbfy2 = np.zeros((len(fx), 1))
    #rbfy3 = np.zeros((len(fx), 1))
    #rbfy4 = np.zeros((len(fx), 1))
    #for i in range(0, len(fx)):
    #    rbfy1[i] = rbfSolution(fx[i], knownParams, coeff1, rbf1, a1)
    #    rbfy2[i] = rbfSolution(fx[i], knownParams, coeff2, rbf2, a2)
    #    rbfy3[i] = rbfSolution(fx[i], knownParams, coeff3, rbf3, a3)
    #    rbfy4[i] = rbfSolution(fx[i], knownParams, coeff4, rbf4, a4)

    """
    fig, ax = plt.subplots()
    rc('text', usetex=True)
    font = {'family' : 'sans-serif',
            'size'   : FONT_SIZE}
    rc('font', **font)
    rc('xtick', labelsize=FONT_SIZE)
    rc('ytick', labelsize=FONT_SIZE)
    ax.plot(fx, fy, 'r-', label=r'$f_{original}$')
    ax.plot(fx, rbfy1, 'b--', label=r'$f_{gaus-RBF}$ mit $a = '+str(a1)+'$')
    ax.plot(fx, rbfy2, 'b:', label=r'$f_{gaus-RBF}$ mit $a = '+str(a2)+'$')
    ax.plot(fx, rbfy3, 'c--', label=r'$f_{imq-RBF}$ mit $a = '+str(a3)+'$')
    ax.plot(fx, rbfy4, 'c:', label=r'$f_{imq-RBF}$ mit $a = '+str(a4)+'$')
    ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)
    ax.legend(loc='lower left', ncol=2)#, mode="expand")
    #ax.legend(loc='lower left', ncol=1)
    ax.set_xlabel('Eingang', fontdict=font)
    ax.set_ylabel('Ausgang', fontdict=font)
    ax.tick_params(labelsize=16., length=6, width=2)
    fig.set_size_inches(8, 5)
    plt.tight_layout()
    plt.savefig('dataOut/radialBasisR2.svg')
    plt.savefig('dataOut/radialBasisR2.pdf')
    plt.show()
    """
