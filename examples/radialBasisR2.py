
import numpy as np

from myLibs.RBF import RBF
from utils.PlotHelper import PlotHelper
from utils.samples import *

if __name__ == '__main__':
    plt1 = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=False)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f_2D,fx))

    plt1.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # now we pretend we only know a view points
    knownParams = [0., 2., 4., 6., 8., 10.]
    knownValues = list(map(f_2D, knownParams))

    fx = fx.reshape((len(fx), 1))
    # RBF 1
    rbf1 = RBF([knownParams], knownValues)
    a1 = 1.
    rbf1.update_param(a1, 'gaus')
    rbfY1 = list(map(rbf1.predict, fx))
    plt1.ax.plot(fx, rbfY1, 'b--', label=r'$\widehat{f}_{gaus-RBF}$ mit $a = ' + str(a1) + '$')

    # RBF 2
    rbf2 = RBF([knownParams], knownValues)
    a2 = 0.24
    rbf2.update_param(a2, 'gaus')
    rbfY2 = list(map(rbf2.predict, fx))
    plt1.ax.plot(fx, rbfY2, 'b:', label=r'$\widehat{f}_{gaus-RBF}$ mit $a = ' + str(a2) + '$')

    # RBF 3
    rbf3 = RBF([knownParams], knownValues)
    a3 = 1.
    rbf3.update_param(a3, 'inverse-multi-quadratic')
    rbfY3 = list(map(rbf3.predict, fx))
    plt1.ax.plot(fx, rbfY3, 'c--', label=r'$\widehat{f}_{imq-RBF}$ mit $a = ' + str(a3) + '$')

    # RBF 4
    rbf4 = RBF([knownParams], knownValues)
    a4 = 1.8
    rbf4.update_param(a4, 'inverse-multi-quadratic')
    rbfY4 = list(map(rbf4.predict, fx))
    plt1.ax.plot(fx, rbfY4, 'c:', label=r'$\widehat{f}_{imq-RBF}$ mit $a = ' + str(a4) + '$')

    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    plt1.finalize(width=6, height=4, legendLoc='lower left', legendNcol=2)
    plt1.ax.legend(loc='upper left', ncol=2)  # , mode="expand")
    plt1.save('../dataOut/radialBasisR2.pdf')
    plt1.show()
