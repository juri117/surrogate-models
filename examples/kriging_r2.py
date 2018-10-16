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
import sys
import os

from mylibs.validation import Validation
from mylibs.structured_sample import StructuredSample

PGF = False

if PGF:
    import matplotlib
    matplotlib.use('pgf')

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/pykriging')
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))+'/../lib/inspyred')
from pyKriging.krige import kriging as PyKriging

#matplotlib.use('pgf')
#pgf_with_custom_preamble = {
#    "pgf.rcfonts": False
#}
#matplotlib.rcParams.update(pgf_with_custom_preamble)

#import matplotlib.pyplot as plt

from mylibs.kriging import Kriging
from myutils.samples import *
from myutils.plot_helper import PlotHelper


class BasinHoppingBoundsLow(object):

    def __init__(self, xmax=[10.], xmin=[0.]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, **kwargs):
        x = kwargs["x_new"]
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin


if __name__ == '__main__':
    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f_2D,fx))

    # now we pretend we only know a view points
    sample = StructuredSample()
    #knownParams = np.array([1., 3., 5., 7., 9., 11.])
    knwonParams = sample.generate_sample_plan(8, 1, [(0., 10.)])
    knownParams = np.array(knwonParams).flatten()
    knownValues = np.array(list(map(f_2D, knownParams)))

    # validate points
    valiParams = np.array([2., 6., 8.])
    valiValues = np.array(list(map(f_2D, valiParams)))
    valiParams = valiParams.reshape((len(valiParams), 1))

    #first fixed exponent here
    p = [1.999966138140631]
    #first fixed factor here
    theta = [0.056389969498335794]

    krig = Kriging(knownParams, knownValues)
    krig2 = PyKriging(knownParams.reshape((len(knownParams), 1)), knownValues)
    krig2.train()

    krig.update_param(theta, p)

    #NegLnLike = calc_likelihood(px, py, theta, p)
    NegLnLike = krig.calc_likelihood()
    print('negLnLike = ' + str(NegLnLike))

    thetas = np.logspace(-2, 3, num=500)
    ps = np.linspace(1., 2., 100)
    likely = np.zeros((len(ps), len(thetas)))
    for it in range(0, len(thetas)):
        for ip in range(0, len(ps)):
            krig.update_param([thetas[it]], [ps[ip]])
            likely[ip][it] = krig.calc_likelihood()
    if True:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.set_xscale('log')
        pcol = ax.pcolor(thetas, ps, likely, cmap='YlOrRd_r')
        fig.colorbar(pcol)
        ax.plot(krig._theta[0], krig._p[0], 'rx')
        ax.set_xlabel('$\theta$')
        ax.set_ylabel('p')

    krig.optimize(opti_algo='grid')
    #krig1.update_param(krig1._theta, krig1._p)

    minLike = krig.calc_likelihood()
    print('minLike = '+str(minLike))
    print('@theta = ' + str(krig._theta[0]))
    print('@p = ' + str(krig._p[0]))

    plt0 = PlotHelper([r'$\theta$', r'Likelihood'], fancy=True, pgf=PGF)
    plt0.ax.semilogx(thetas, likely[-1])
    plt0.ax.semilogx(krig._theta[0], minLike, 'r+', markersize=8, label='Minimum')
    plt0.finalize(width=4, height=2, legendLoc='upper right', legendNcol=1)
    plt0.save('../data_out/plot/krigingR2likelihood.pdf')
    #plt0.show()

    ###### validate
    vali = Validation()
    vali_r = vali.run_full_analysis(fx.reshape((len(fx), 1)), fy,
                                    knownParams.reshape((len(knownParams), 1)), knownValues,
                                    valiParams, valiValues,
                                    krig.predict, Kriging)
    print('avg deviation: {:.3e} (-> {:.3f}%)'.format(vali_r.deviation, vali_r.deviation * 100.))
    print('rmse: {:f}'.format(vali_r.rmse))
    print('mae: {:f}'.format(vali_r.mae))
    print('rae: {:s}'.format(str(vali_r.rae)))
    print('press: {:f}'.format(vali_r.press))

    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=PGF)
    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')
    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    krigY = list(map(krig.predict, fx.reshape((len(fx), 1))))
    plt1.ax.plot(fx, krigY, 'b--', label=r'$\widehat{f}_{krig}$ mit $\theta = ' +'{0:.3f}'.format(krig._theta[0]) + '$, $p = ' + '{0:.1f}'.format(krig._p[0]) + '$')

    # plot PyKriging
    #krig2Y = list(map(krig2.predict, fx.reshape((len(fx), 1))))
    #plt1.ax.plot(fx, krig2Y, '--', color='green', label=r'$\widehat{f}_{PyKrig}$')

    # scipy minimize
    res = minimize(krig.predict, [3.], method='SLSQP', bounds=[(knownParams[0], knownParams[-1])])
    plt1.ax.plot([res.x], [krig.predict(res.x)], 'co', label=r'Minimum SLSQP')
    print('SLSQP tries: {:d}'.format(res.nit))

    res2 = optimize.differential_evolution(krig.predict, [(knownParams[0], knownParams[-1])], disp=True)
    plt1.ax.plot([res2.x], [krig.predict(res2.x)], 'go', label=r'Minimum, diff. evo.')
    print('differential_evolution tries: {:d}'.format(res2.nit))

    #bounds = BasinHoppingBoundsLow(xmax=[knownParams[-1]], xmin=[knownParams[0]])
    #minimizer_kwargs = dict(method='SLSQP', bounds=[(knownParams[0], knownParams[-1])], options={'disp': False})
    #res3 = optimize.basinhopping(krig.predict, [3.], disp=True, minimizer_kwargs=minimizer_kwargs, accept_test=bounds)
    #plt1.ax.plot([res3.x], [krig.predict(res3.x)], 'go', label=r'Minimum, Basin-hop.')
    #print('differential_evolution tries: {:d}'.format(res3.nit))

    plt1.finalize(width=6, height=4, legendLoc='upper left', legendNcol=1)
    plt1.save('../data_out/plot/krigingR2.pdf')
    plt1.show()
