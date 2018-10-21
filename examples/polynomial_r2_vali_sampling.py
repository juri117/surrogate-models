
import numpy as np

from mylibs.polynomial import Polynomial
from myutils.plot_helper import PlotHelper
from myutils.samples import *
from mylibs.validation import Validation
from mylibs.validation import ValidationResults
from mylibs.latin_hyper_cube import LatinHyperCube
from mylibs.halton import Halton
from mylibs.structured_sample import StructuredSample

if __name__ == '__main__':
    polyPlot = PlotHelper(['Eingang', 'Ausgang'], fancy=True, pgf=True)

    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f_2D,fx))

    polyPlot.ax.plot(fx, fy, 'r-', label=r'$f_{original}$')

    # validate points
    valiParams = np.array([1., 5. , 9.])
    valiValues = np.array(list(map(f_2D, valiParams)))
    valiParams = valiParams.reshape((len(valiParams), 1))

    fx = fx.reshape((len(fx), 1))

    vali_res = []
    poly_plot_points = [2, 3, 4, 5, 6, 7]
    sampling_points = [2, 3, 4, 5, 6, 7, 8, 9, 10] # [2, 3, 4, 5, 6, 7]
    plot_styles = ['k', 'k', 'k', 'silver', 'm', 'darksalmon', 'limegreen', 'cornflowerblue', 'b', 'm', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b', 'b']
    for n in sampling_points:

        # now we pretend we only know a view points
        sample = StructuredSample()
        knownParams = sample.generate_sample_plan(n, 1, [(0., 10.)])
        knownParams = np.array(knownParams).flatten()
        #knownParams = np.array([0., 2., 4., 6., 8., 10.])
        knownValues = np.array(list(map(f_2D, knownParams)))

        # Polynomials
        poly = Polynomial(knownParams, knownValues)

        print('##############')
        print('Order: {:d}'.format(n))
        o = 6
        poly.update_param(o)
        rbfY1 = list(map(poly.predict, fx))
        if n in poly_plot_points:
            polyPlot.ax.scatter(knownParams, knownValues, facecolors='none', edgecolors=plot_styles[n], s=100)
            polyPlot.ax.plot(fx, rbfY1, '--', color=plot_styles[n], label=r'$\widehat{f}_{polyn}$, $n = ' + str(n) + '$')
        vali = Validation()

        vali = Validation()
        vali_r = vali.run_full_analysis(fx, fy,
                                        knownParams.reshape((len(knownParams), 1)), knownValues,
                                        valiParams, valiValues,
                                        poly.predict, Polynomial,
                                        update_params=[o])
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

    #polyPlot.ax.plot(knownParams, knownValues, 'ro', label=r'St\"utzstellen', markersize=10)

    polyPlot.ax.scatter(valiParams, valiValues, facecolors='none', edgecolors='r', label=r'Vali.-Punkte', s=100)

    polyPlot.finalize(width=3.9, height=4, legendLoc='upper center', legendNcol=2, bbox_to_anchor=(0.5, -0.35))
    polyPlot.ax.set_ylim([2, 10])
    import matplotlib.pyplot as plt
    plt.subplots_adjust(bottom=0.494)
    polyPlot.save('../data_out/plot/polynomialR2sampling.pdf')


    # plot validation
    valiPlot = PlotHelper(['Stützstellen', 'Wert'], fancy=True, pgf=True)
    valiPlot.ax.plot(sampling_points, [r.deviation for r in vali_res], '--', color='silver', label=r'$\O$ -Abweichung in $\%$')
    valiPlot.ax.plot(sampling_points, [r.rmse for r in vali_res], '-', color='orange', label=r'$RMSE$')
    valiPlot.ax.plot(sampling_points, [r.mae for r in vali_res], '-', color='green', label=r'$MAE$')
    valiPlot.ax.plot(sampling_points, np.array([r.press for r in vali_res]) * 0.1, 'k-', label=r'$PRESS \cdot 0.1$')
    valiPlot.finalize(width=2.8, height=4, legendLoc='upper center', bbox_to_anchor=(0.5, -0.25))
    valiPlot.ax.set_xlim([2, 10])
    valiPlot.ax.set_ylim([0, 2.5])
    plt.locator_params(nbins=4)
    plt.subplots_adjust(bottom=0.42)
    plt.subplots_adjust(left=0.22)
    polyPlot.save('../data_out/plot/polynomialR2valiSampling.pdf')
    valiPlot.show()


