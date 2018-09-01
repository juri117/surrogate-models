__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :kriging in R2 with simple function
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
from scipy.optimize import minimize
from scipy import optimize
import matplotlib

PGF = False

#matplotlib.use('pgf')
#pgf_with_custom_preamble = {
#    "pgf.rcfonts": False
#}
#matplotlib.rcParams.update(pgf_with_custom_preamble)

#import matplotlib.pyplot as plt

from myLibs.Kriging import Kriging
from utils.samples import *
from utils.PlotHelper import PlotHelper


if __name__ == '__main__':
    # the smooth whole function
    fx = np.linspace(1, 11, 1001)
    fy = list(map(f_2D,fx))

    # now we pretend we only know a view points
    px = [1., 3., 5., 7., 9., 11.]
    py = list(map(f_2D,px))
    #first fixed exponent here
    p = [2.]
    #first fixed factor here
    theta = [.5]

    krig1 = Kriging(px, py)
    krig1.update_param(theta, p)

    #NegLnLike = calc_likelihood(px, py, theta, p)
    NegLnLike = krig1.calc_likelihood()
    print('negLnLike = ' + str(NegLnLike))

    #thetas = np.linspace(0.01, 10, 1000+1)
    thetas = np.logspace(-2, 3, num=500)
    ps = np.linspace(1., 2., 100)
    likely = np.zeros((len(ps), len(thetas)))
    for it in range(0, len(thetas)):
        for ip in range(0, len(ps)):
            krig1.update_param([thetas[it]], [ps[ip]])
            likely[ip][it] = krig1.calc_likelihood()

    krig1.optimize()
    #krig1.update_param(krig1._theta, krig1._p)

    minLike = krig1.calc_likelihood()
    print('minLike = '+str(minLike))
    print('@theta = ' + str(krig1._theta[0]))
    print('@p = ' + str(krig1._p[0]))

    plt0 = PlotHelper([r'$\theta$', r'Likelihood'], fancy=False, pgf=PGF)
    plt0.ax.semilogx(thetas, likely[-1])
    plt0.ax.semilogx(krig1._theta[0], minLike, 'rx', markersize=10, label='Minimum')
    plt0.finalize(width=6, height=3.5, legendLoc='upper right', legendNcol=1)
    plt0.save('../dataOut/krigingR2likelihood.pdf')
    #plt0.show()

    #fig, ax = plt.subplots()
    #ax.set_xscale('log')
    #pcol = ax.pcolor(thetas, ps, likely, cmap='viridis_r')
    #fig.colorbar(pcol)
    #ax.plot(krig1._theta[0], krig1._p[0], 'rx')
    #ax.set_xlabel('$\theta$')
    #ax.set_ylabel('p')
    #plt.show()

    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=PGF)

    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')
    plt1.ax.plot(px, py, 'ro', label=r'St\"utzstellen', markersize=10)

    krigY = list(map(krig1.predict, fx.reshape((len(fx), 1))))
    plt1.ax.plot(fx, krigY, 'b-', label=r'$\widehat{f}_{krig}$ mit $\theta = '+'{0:.3f}'.format(krig1._theta[0])+'$, $p = '+'{0:.1f}'.format(krig1._p[0])+'$')

    # scipy minimize
    res = minimize(krig1.predict, [3.], method='SLSQP', bounds=[(px[0], px[-1])])
    plt1.ax.plot([res.x], [krig1.predict(res.x)], 'co', label=r'Minimum SLSQP')

    res = optimize.differential_evolution(krig1.predict, [(px[0], px[-1])])
    plt1.ax.plot([res.x], [krig1.predict(res.x)], 'go', label=r'Minimum, diff. evo.')

    plt1.finalize(width=6, height=4, legendLoc='upper left', legendNcol=1)
    plt1.save('../dataOut/krigingR2.pdf')
    plt1.show()
