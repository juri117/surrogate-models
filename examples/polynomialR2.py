
import numpy as np

from myLibs.Polynomial import Polynomial
from utils.PlotHelper import PlotHelper
from utils.samples import *
from myLibs.Validation import Validation
from myLibs.Validation import ValidationResults
from myLibs.LatinHyperCube import LatinHyperCube
from myLibs.Halton import Halton
from myLibs.StructuredSample import StructuredSample

if __name__ == '__main__':
    polyPlot = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=False)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f_2D,fx))

    polyPlot.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # now we pretend we only know a view points
    sample = StructuredSample()
    knwonParams = sample.generate_sample_plan(6, 1, [(0., 10.)])
    knownParams = np.array(knwonParams).flatten()
    #knownParams = np.array([0., 2., 4., 6., 8., 10.])
    knownValues = np.array(list(map(f_2D, knownParams)))

    # validate points
    valiParams = np.array([1., 5. , 9.])
    valiValues = np.array(list(map(f_2D, valiParams)))
    valiParams = valiParams.reshape((len(valiParams), 1))

    fx = fx.reshape((len(fx), 1))
    # Polynomials
    poly = Polynomial(knownParams, knownValues)
    vali_res = []
    poly_plot_points = [1, 3, 5, 6, 7]
    degrees = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10] #[3, 5, 6, 7]
    plot_styles = ['k', 'silver', 'k', 'darksalmon', 'k', 'cornflowerblue', 'limegreen', 'm', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k', 'k']
    for i in degrees:
        print('##############')
        print('Order: {:d}'.format(i))
        o = i
        poly.update_param(o)
        rbfY1 = list(map(poly.predict, fx))
        if i in poly_plot_points:
            polyPlot.ax.plot(fx, rbfY1, '--', color=plot_styles[i], label=r'$\widehat{f}_{polyn}$ mit $o = ' + str(o) + '$')

        vali = Validation()
        vali_r = vali.run_full_analysis(fx, fy,
                                        knownParams.reshape((len(knownParams), 1)), knownValues,
                                        valiParams, valiValues,
                                        poly.predict, Polynomial,
                                        update_params=[i])
        vali_res.append(vali_r)
        #deviation = vali.calc_deviation(fx, fy, poly.predict)
        #rmse = vali.calc_rmse(valiParams, valiValues, poly.predict)
        #mae = vali.calc_mae(valiParams, valiValues, poly.predict)
        #rae = vali.calc_rae(valiParams, valiValues, poly.predict)
        #press = vali.calc_press(knownParams.reshape((len(knownParams), 1)), knownValues, poly.predict, Polynomial, update_params=[i])

        print('avg deviation: {:.3e} (-> {:.3f}%)'.format(vali_r.deviation, vali_r.deviation * 100.))
        print('rmse: {:f}'.format(vali_r.rmse))
        print('mae: {:f}'.format(vali_r.mae))
        print('rae: {:s}'.format(str(vali_r.rae)))
        print('press: {:f}'.format(vali_r.press))

    polyPlot.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)
    polyPlot.ax.scatter(valiParams, valiValues, facecolors='none', edgecolors='r', label=r'Validierungspunkte', s=100)
    polyPlot.finalize(width=6, height=4, legendLoc='upper left', legendNcol=2)
    polyPlot.save('../dataOut/plot/polynomialR2.pdf')


    # plot validation
    valiPlot = PlotHelper(['Grad', 'Wert'], fancy=True, pgf=True)
    valiPlot.ax.plot(degrees, [r.deviation for r in vali_res], '--', color='silver', label=r'\O Abweichung in $\%$')
    valiPlot.ax.plot(degrees, [r.rmse for r in vali_res], '-', color='orange', label=r'$RMSE$')
    valiPlot.ax.plot(degrees, [r.mae for r in vali_res], '-', color='green', label=r'$MAE$')
    valiPlot.ax.plot(degrees, np.array([r.press for r in vali_res]) * 0.1, 'k-', label=r'$PRESS * 0.1$')
    valiPlot.finalize(width=4.5, height=3, legendLoc='upper left')
    valiPlot.ax.set_ylim([0, 2.5])
    polyPlot.save('../dataOut/plot/polynomialR2vali.pdf')
    valiPlot.show()


