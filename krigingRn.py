
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize
import time

from Kriging import Kriging
from KrigingForrester import KrigingForrester
from utils.PlotHelper import PlotHelper
from utils.TimeTrack import TimeTrack

FONT_SIZE=14
font = {'family':'sans-serif', 'size':FONT_SIZE}

def printMat(mat, octave=False):
    if octave:
        print('['+';\n'.join([''.join(['{:10.3} '.format(item) for item in row]) for row in mat])+']')
    else:
        print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x, y):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4 + 0.05*y**2 - 0.001*y**4 - 0.005*y**3


if __name__ == '__main__':
    t1 = TimeTrack('OverAllTimer')
    # the smooth whole function
    fx = np.linspace(-2, 12, 101)
    fy = np.linspace(-2, 12, 101)
    fz = np.zeros((len(fx), len(fy)))
    for iX in range(0,len(fx)):
        for iY in range(0,len(fy)):
            fz[iX][iY] = f(fx[iX], fy[iY])

    # now we pretend we only know a view points
    knownParams = []
    pxEdge = [0., 2., 4., 6., 8., 10.]
    pyEdge = [0., 2., 4., 6., 8., 10.]
    px = []
    py = []
    pz = []
    for iX in range(0, len(pxEdge)):
        for iY in range(0, len(pyEdge)):
            px.append(pxEdge[iX])
            py.append(pyEdge[iY])
            pz.append(f(pxEdge[iX], pyEdge[iY]))
    knownParams.append(px)
    knownParams.append(py)
    knownParams = np.array(knownParams)
    knownValues = np.array(pz)

    krig = Kriging(knownParams, knownValues)

    p = [1.8, 1.8]

    krig.update_param([0.001, 0.001], p)
    print(str(krig.calc_likelihood()))
    #print(str(krig.calc_likelihood_v2([0.001, 0.001], p)))

    #thetas = np.linspace(0.001, 0.1, 100 + 1)
    thetas = np.logspace(-3, 1, num=20)
    likelyX1 = []
    likelyX2 = []
    likely = np.zeros((len(thetas), len(thetas)))
    for i1 in range(0, len(thetas)):
        #print(str(i1) + ' / ' + str(len(thetas)))
        for i2 in range(0, len(thetas)):
            krig.update_param([thetas[i1], thetas[i2]], p)
            likely[i2][i1] = krig.calc_likelihood()

    krig.optimize()
    minLike = krig.calc_likelihood()
    print('minLike = ' + str(minLike))
    print('@theta1 = ' + str(krig._theta[0]))
    print('@theta2 = ' + str(krig._theta[1]))

    pltTheta = PlotHelper([r'$\theta_{1}$', r'$\theta_{1}$'], fancy=False)
    pltTheta.ax.set_xscale('log')
    pltTheta.ax.set_yscale('log')
    pcol = pltTheta.ax.pcolor(thetas, thetas, likely)
    pltTheta.fig.colorbar(pcol)
    pltTheta.ax.plot(krig._theta[0], krig._theta[1], 'rx')
    pltTheta.finalize()
    pltTheta.show()

    #plot original function and points
    plotX, plotY = np.meshgrid(fx, fy)
    plt1 = PlotHelper(['Eingang 1', 'Eingang 2', 'Ausgang'], fancy=False)

    surf1 = plt1.ax.plot_wireframe(plotX, plotY, fz, color='r', label=r'$f_{original}$', rcount=20, ccount=20, linewidths=1,
                              alpha=0.5)  # , rstride=1, cstride=1)#, cmap=cm.coolwarm) # ,linewidth=0, antialiased=False

    scat1 = plt1.ax.scatter(py, px, pz, c='r', marker='o', s=10, label=r'St\"utzstellen')



    krigingSol = np.zeros((len(fx), len(fy)))
    for iX in range(0, len(fx)):
        #print(str(iX) + ' of ' + str(len(fx)))
        for iY in range(0, len(fy)):
            coords = [fx[iX], fy[iY]]
            krigingSol[iX][iY] = krig.predict(coords)

    surf2 = plt1.ax.plot_wireframe(plotX, plotY, krigingSol, color='b', label=r'$f_{kriging}$', rcount=20, ccount=20, linewidths=1,
                              alpha=0.5)#, antialiased=True)

    plt1.ax.view_init(20, 50)
    plt1.finalize()
    #plt1.save('dataOut/krigingRn.svg')
    #plt1.save('dataOut/krigingRn.pdf')
    #plt1.animate()
    t1.toc()
    plt1.show()
    print('done')














