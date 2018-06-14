
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import scipy

def printMat(mat):
    print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4

# gaus RBF (GS)
def gausRBF(a, r):
    #return r
    #return math.e**(-((a*r)**2))
    return math.sqrt(1 + (a*r)**2)

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

# the smooth whole function
fx = np.linspace(0, 10, 1001)
fy = list(map(f,fx))

# now we pretend we only know a view points
px = [0., 2., 4., 6., 8., 10.]
py = list(map(f,px))

a1 = 1.
coeff1 = rbf_calc_coefficiants(a1, px, py, gausRBF)

a2 = 0.17
coeff2 = rbf_calc_coefficiants(a2, px, py, gausRBF)

rbfy1 = np.zeros((len(fx), 1))
rbfy2 = np.zeros((len(fx), 1))
for i in range(0, len(fx)):
    rbfy1[i] = rbfSolution(fx[i], px, coeff1, gausRBF, a1)
    rbfy2[i] = rbfSolution(fx[i], px, coeff2, gausRBF, a2)

#rc('text', usetex=True)
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.plot(fx, fy, 'b-', label=r'f_{original}')
plt.plot(fx, rbfy1, 'g-', label=r'f_{rbf} mit a = 1')
plt.plot(fx, rbfy2, 'r--', label=r'f_{rbf} mit a = 0.17')
plt.plot(px, py, 'r*', label=r'Stuetzstellen')
plt.legend()

plt.show()

