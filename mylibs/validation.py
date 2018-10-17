__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :validation of surrogates
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
import math
import sys

from myutils.plot_helper import PlotHelper

class Validation:

    def __init__(self):
        pass

    '''
    calc the deviation of a matrix with known solutions to the surrogate solution
    '''
    def calc_deviation(self, params, values, surro_func):
        count = 0
        sum_deviation = 0
        # sample_indices = np.array([known_x_i, known_y_i]).T.tolist()
        for i in range(0, len(params)):
            #for i_s in range(0, len(y)):
            devi = values[i] - surro_func(params[i])
            sum_deviation += abs(devi)
            count += 1
        avg_deviation = sum_deviation / count
        avg_deviation_per = avg_deviation / np.array(values).mean()
        return avg_deviation_per

    '''
    :return root mean square error (RMSE)
    '''
    def calc_rmse(self, vali_x, vali_fx, surro_func):
        sum = 0.
        for i in range(0, len(vali_x)):
            res = surro_func(vali_x[i])
            sum += (vali_fx[i] - res) ** 2
        return math.sqrt(sum / len(vali_x))

    '''
    :return maximum absolute error (MAE)
    '''
    def calc_mae(self, vali_x, vali_fx, surro_func):
        mae = -1 * float('inf')
        for i in range(0, len(vali_x)):
            res = surro_func(vali_x[i])
            ae = abs(vali_fx[i] - res)
            mae = max(mae, ae)
        return mae

    '''
    :return relative absolute error (RAE)
    '''
    def calc_rae(self, vali_x, vali_fx, surro_func):
        pred = list(map(surro_func, vali_x))
        return abs(np.divide((vali_fx - pred), vali_fx))

    '''
    :return prediction sum of squares (PRESS)
    '''
    def calc_press(self, known_x, known_fx, surro_func, surro_class, update_params=None):
        sum = 0.
        for i in range(0, len(known_x)):
            res = surro_func(known_x[i])
            x = np.array(known_x.copy())
            x = np.delete(x, 0, axis=0)
            fx = np.array(known_fx.copy())
            fx = np.delete(fx, 0, axis=0)
            try:
                sur = surro_class(x, fx)
                if update_params != None:
                    if len(update_params) == 1:
                        sur.update_param(update_params[0])
                    elif len(update_params) == 2:
                        sur.update_param(update_params[0], update_params[1])
                    else:
                        print('Validation.calc_press please update to support more params than 2')
                sur.train()
            except Exception as e:
                print('WARNING: missing train() method: ' + str(e))
                return 0.
            res2 = sur.predict(known_x[i])
            sum += (res - res2)**2.
        return sum / len(known_x)

    def run_full_analysis(self, params, values, known_x, known_fx, vali_x, vali_fx, surro_func, surro_class, update_params=None):
        res = ValidationResults()
        res.deviation = self.calc_deviation(params, values, surro_func)
        res.rmse = self.calc_rmse(vali_x, vali_fx, surro_func)
        res.mae = self.calc_mae(vali_x, vali_fx, surro_func)
        res.rae = self.calc_rae(vali_x, vali_fx, surro_func)
        res.press = self.calc_press(known_x, known_fx, surro_func, surro_class, update_params=update_params)
        return res

    '''
    def plot_derivation2d(self, xs, ys, vals, surro_func, sample_x=None, sample_y=None, opti_x=None, opti_y=None):
        deri_plot = PlotHelper(['param1', 'param2'], fancy=False, pgf=False)
        dev = np.zeros(vals.shape)
        for xi in range(0, len(xs)):
            for yi in range(0, len(ys)):
                devi = (abs(vals[yi][xi] - surro_func([xs[xi], ys[yi]])) / np.array(vals).mean()) *100.
                dev[yi][xi] = devi
        pcol = deri_plot.ax.pcolor(xs, ys, dev, cmap='YlOrRd')
        pcol.set_clim(0, 5.)
        cbar = deri_plot.fig.colorbar(pcol)
        show_legend = False
        if sample_x.all() != None and sample_y.all() != None:
            deri_plot.ax.plot(sample_x, sample_y, 'bo', label='sampling points')
            show_legend = True
        deri_plot.finalize(width=6, height=5, legendLoc='upper right', show_legend=show_legend)
    '''


class ValidationResults():
    def __init__(self):
        self.deviation = 0.
        self.rmse = 0.
        self.mae = 0.
        self.rae = 0.
        self.press = 0.


if __name__ == '__main__':
    val = Validation()


