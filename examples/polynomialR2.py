
import numpy as np

from myLibs.Polynomial import Polynomial
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
    # Polynomials
    poly = Polynomial([knownParams], knownValues)

    for i in [3, 5, 6, 7]:
        o = i
        poly.update_param(o)
        rbfY1 = list(map(poly.predict, fx))
        plt1.ax.plot(fx, rbfY1, '--', label=r'$\widehat{f}_{polyn}$ mit $o = ' + str(o) + '$')

    plt1.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    plt1.finalize(width=6, height=4, legendLoc='lower left', legendNcol=2)
    plt1.ax.legend(loc='upper left', ncol=2)  # , mode="expand")
    plt1.save('../dataOut/plot/polynomialR2.pdf')
    plt1.show()
