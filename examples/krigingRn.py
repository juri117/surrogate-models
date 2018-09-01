__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :kriging in Rn
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np

from myLibs.Kriging import Kriging
from utils.PlotHelper import PlotHelper
from utils.TimeTrack import TimeTrack
from utils.samples import *


if __name__ == '__main__':
    t1 = TimeTrack('OverAllTimer')
    plt1 = PlotHelper(['Eingang 1', 'Eingang 2', 'Ausgang'], fancy=False)

    # now we pretend we only know a view points
    fx = np.linspace(-2, 12, 101)
    fy = np.linspace(-2, 12, 101)
    plt1.plot_function_3D(f_3D, fx, fy, r'$f_{original}$', color='r')
    # the smooth whole function

    # now we pretend we only know a view points
    pxEdge = [0., 2., 4., 6., 8., 10.]
    pyEdge = [0., 2., 4., 6., 8., 10.]
    px, py, knownValues = generate_sample_data(f_3D, pxEdge, pyEdge)
    knownParams = []
    knownParams.append(px)
    knownParams.append(py)

    scat1 = plt1.ax.scatter(px, py, knownValues, c='r', marker='o', s=10, label=r'St\"utzstellen')

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
    #pltTheta.show()

    plt1.plot_function_3D(krig.predict, fx, fy, r'$\widehat{f}_{krig}$', color='b')
    plt1.ax.view_init(20, 50)
    plt1.finalize()
    #plt1.save('dataOut/krigingRn.svg')
    #plt1.save('dataOut/krigingRn.pdf')
    #plt1.animate()
    t1.toc()
    plt1.show()
    print('done')
