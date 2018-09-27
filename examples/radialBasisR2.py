
import numpy as np

from myLibs.RBF import RBF
from myLibs.StructuredSample import StructuredSample
from myLibs.Validation import Validation
from myUtils.PlotHelper import PlotHelper
from myUtils.samples import *

if __name__ == '__main__':
    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=False, pgf=False)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = np.array(list(map(f_2D,fx)))

    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # now we pretend we only know a view points
    sample = StructuredSample()
    knwonParams = sample.generate_sample_plan(6, 1, [(1., 11.)])
    knownParams = np.array(knwonParams).flatten()
    knownValues = np.array(list(map(f_2D, knownParams)))

    fx = fx.reshape((len(fx), 1))

    aas = [1., 0.24, 1., 1.8]
    rbfs = ['gaus', 'gaus', 'imq', 'imq']
    plot_styles = ['b--', 'b:', 'c--', 'c:']

    for i in range(0, len(aas)):
        print('##############')
        print('RBF({:d})'.format(i))
        rbf = RBF(knownParams, knownValues)
        rbf.update_param(aas[i], rbfs[i])
        rbfY1 = list(map(rbf.predict, fx))
        plt1.ax.plot(fx, rbfY1, plot_styles[i], label=r'$\widehat{f}_{' + rbfs[i] + '-RBF}$ mit $a = ' + str(aas[i]) + '$')

        vali = Validation()
        deviation = vali.calc_deviation(fx.reshape((len(fx), 1)), fy, rbf.predict)
        rmse = vali.calc_rmse(knownParams.reshape((len(knownParams), 1)), knownValues, rbf.predict)
        mae = vali.calc_mae(knownParams.reshape((len(knownParams), 1)), knownValues, rbf.predict)
        rae = vali.calc_rae(knownParams.reshape((len(knownParams), 1)), knownValues, rbf.predict)
        press = vali.calc_press(knownParams.reshape((len(knownParams), 1)), knownValues, rbf.predict, RBF, update_params=[aas[i], rbfs[i]])
        print('avg deviation: {:.3e} (-> {:.3f}%)'.format(deviation, deviation * 100.))
        print('rmse: {:f}'.format(rmse))
        print('mae: {:f}'.format(mae))
        print('rae: {:s}'.format(str(rae)))
        print('press: {:f}'.format(press))

    '''
    # RBF 1
    rbf1 = RBF(knownParams, knownValues)
    a1 = 1.
    rbf1.update_param(a1, 'gaus')
    rbfY1 = list(map(rbf1.predict, fx))
    plt1.ax.plot(fx, rbfY1, 'b--', label=r'$\widehat{f}_{gaus-RBF}$ mit $a = ' + str(a1) + '$')

    # RBF 2
    rbf2 = RBF(knownParams, knownValues)
    a2 = 0.24
    rbf2.update_param(a2, 'gaus')
    rbfY2 = list(map(rbf2.predict, fx))
    plt1.ax.plot(fx, rbfY2, 'b:', label=r'$\widehat{f}_{gaus-RBF}$ mit $a = ' + str(a2) + '$')

    # RBF 3
    rbf3 = RBF(knownParams, knownValues)
    a3 = 1.
    rbf3.update_param(a3, 'inverse-multi-quadratic')
    rbfY3 = list(map(rbf3.predict, fx))
    plt1.ax.plot(fx, rbfY3, 'c--', label=r'$\widehat{f}_{imq-RBF}$ mit $a = ' + str(a3) + '$')

    # RBF 4
    rbf4 = RBF(knownParams, knownValues)
    a4 = 1.8
    rbf4.update_param(a4, 'inverse-multi-quadratic')
    rbfY4 = list(map(rbf4.predict, fx))
    plt1.ax.plot(fx, rbfY4, 'c:', label=r'$\widehat{f}_{imq-RBF}$ mit $a = ' + str(a4) + '$')
    '''

    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    plt1.finalize(width=6, height=4, legendLoc='lower left', legendNcol=2)
    plt1.ax.legend(loc='upper left', ncol=2)  # , mode="expand")
    plt1.save('../dataOut/plot/radialBasisR2.pdf')
    plt1.show()
