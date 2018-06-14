
import numpy as np
import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from matplotlib import rc
import scipy

def printMat(mat):
    print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x, y):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4 + 0.05*y**2 - 0.001*y**4 - 0.005*y**3

# gaus RBF (GS)
def gausRBF(a, r):
    #return r
    #return math.e**(-((a*r)**2))
    return math.sqrt(1 + (a*r)**2)

# todo: other RBFunctions here...

def rbf_calc_coefficiants(rbf_const, knownCoord, knownVal, rbf):
    paramCount, knownCount = np.array(knownCoord).shape
    mat = np.zeros((knownCount, knownCount))
    radialMat = np.zeros((knownCount, knownCount))
    for iRow in range(0, knownCount):
        for iColumn in range(iRow, knownCount):
            radSum = 0.
            for iParam in range(0,paramCount):
                radSum += (knownCoord[iParam][iRow] - knownCoord[iParam][iColumn])**2.
            radius = math.sqrt(radSum)
            radialMat[iRow][iColumn] = radius
            mat[iRow][iColumn] = rbf(rbf_const, radius)
            mat[iColumn][iRow] = mat[iRow][iColumn]
    #print('radialMat:')
    #printMat(radialMat)
    print('mat:')
    printMat(mat)
    return np.linalg.solve(mat, knownVal)

def rbfSolution(x, knownCoord, coefficients, rbf, rbf_const):
    paramCount, knownCount = np.array(knownCoord).shape
    n = len(knownCoord)
    res = 0.
    for i in range(0, knownCount):
        radSum = 0.
        for iParam in range(0, paramCount):
            radSum += (x[iParam] - knownCoord[iParam][i]) ** 2.
        radius = math.sqrt(radSum)
        res += coefficients[i] * rbf(rbf_const, radius)
    return res


fig = plt.figure()
ax = fig.gca(projection='3d')

# the smooth whole function
fx = np.linspace(-2, 12, 201)
fy = np.linspace(-2, 12, 201)
fz = np.zeros((len(fx), len(fy)))
for iX in range(0,len(fx)):
    for iY in range(0,len(fy)):
        fz[iX][iY] = f(fx[iX], fy[iY])

plotX, plotY = np.meshgrid(fx, fy)
fig = plt.figure()
ax = fig.gca(projection='3d')
surf1 = ax.plot_wireframe(plotX, plotY, fz, rcount=20, ccount=20)#, rstride=1, cstride=1)#, cmap=cm.coolwarm) # ,linewidth=0, antialiased=False
#plt.show()





# now we pretend we only know a view points
knownParams = []
pxEdge = [0., 2., 4., 6., 8., 10.]
pyEdge = [0., 2., 4., 6., 8., 10.]
px=[]
py=[]
pz=[]
for iX in range(0,len(pxEdge)):
    for iY in range(0,len(pyEdge)):
        px.append(pxEdge[iX])
        py.append(pyEdge[iY])
        pz.append(f(pxEdge[iX],pyEdge[iY]))
knownParams.append(px)
knownParams.append(py)

ax.scatter(py, px, pz, c='g', marker='o')

a1 = 1.
coeff1 = rbf_calc_coefficiants(a1, knownParams, pz, gausRBF)

a2 = 0.17
coeff2 = rbf_calc_coefficiants(a2, knownParams, pz, gausRBF)

rbfz1 = np.zeros((len(fx), len(fy)))
rbfz2 = np.zeros((len(fx), len(fy)))
for iX in range(0,len(fx)):
    for iY in range(0,len(fy)):
        coords = [fx[iX], fy[iY]]
        rbfz1[iX][iY] = rbfSolution(coords, knownParams, coeff1, gausRBF, a1)
        rbfz2[iX][iY] = rbfSolution(coords, knownParams, coeff2, gausRBF, a2)

surf2 = ax.plot_wireframe(plotX, plotY, rbfz1, color='r', rcount=20, ccount=20)#, rstride=1, cstride=1)#, cmap=cm.coolwarm) # ,linewidth=0, antialiased=False
plt.show()


