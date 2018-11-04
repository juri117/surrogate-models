
import numpy as np

from mylibs.rbf import RBF
from mylibs.structured_sample import StructuredSample
from mylibs.validation import Validation
from myutils.plot_helper import PlotHelper
from myutils.samples import *

if __name__ == '__main__':
    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=False)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = np.array(list(map(f_2d, fx)))

    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # now we pretend we only know a view points
    sample = StructuredSample()
    knwonParams = sample.generate_sample_plan(6, 1, [(0., 10.)])
    knownParams = np.array(knwonParams).flatten()
    knownValues = np.array(list(map(f_2d, knownParams)))

    fx = fx.reshape((len(fx), 1))

    aas = [1., 0.24, 1., 0.24]
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

    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    plt1.finalize(width=6, height=4, legendLoc='lower left', legendNcol=2, tighten_layout=True)
    plt1.ax.legend(loc='upper left', ncol=2)  # , mode="expand")
    plt1.save('../data_out/plot/radialBasisR2.pdf')
    plt1.show()
