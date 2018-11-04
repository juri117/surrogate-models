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


class Validation:

    def __init__(self):
        pass

    def calc_deviation(self, params, values, surro_func):
        """
        :param list of entries for fem grid calculation
        :param values: list of results of params
        :param surro_func: pointer to the surrogates predict function
        :return: the deviation of a matrix with known solutions to the surrogate solution
        """
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

    def calc_rmse(self, vali_x, vali_fx, surro_func):
        """
        :param vali_x: list of list of entry values
        :param vali_fx: list of results to each vali_x
        :param surro_func: pointer to the surrogates predict function
        :return: root mean square error (RMSE)
        """
        sum = 0.
        for i in range(0, len(vali_x)):
            res = surro_func(vali_x[i])
            sum += (vali_fx[i] - res) ** 2
        return math.sqrt(sum / len(vali_x))

    def calc_mae(self, vali_x, vali_fx, surro_func):
        """
        :param vali_x: list of list of entry values
        :param vali_fx: list of results to each vali_x
        :param surro_func: pointer to the surrogates predict function
        :return: maximum absolute error (MAE)
        """
        mae = -1 * float('inf')
        for i in range(0, len(vali_x)):
            res = surro_func(vali_x[i])
            ae = abs(vali_fx[i] - res)
            mae = max(mae, ae)
        return mae

    def calc_rae(self, vali_x, vali_fx, surro_func):
        """
        :param vali_x: list of list of entry values
        :param vali_fx: list of results to each vali_x
        :param surro_func: pointer to the surrogates predict function
        :return: list of relative absolute error (RAE) for each sample point
        """
        pred = list(map(surro_func, vali_x))
        return abs(np.divide((vali_fx - pred), vali_fx))

    def calc_press(self, known_x, known_fx, surro_func, surro_class, update_params=None):
        """
        :param known_x: list of sampling points (matrix)
        :param known_fx: list of results for known_x
        :param surro_func: pointer to the surrogates predict function
        :param surro_class: class of the used surrogate
        :param update_params: parameters to pass to surrogate on fitting
        :return: prediction sum of squares (PRESS)
        """
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
        return math.sqrt(sum / len(known_x))

    def run_full_analysis(self, params, values, known_x, known_fx, vali_x, vali_fx, surro_func, surro_class, update_params=None):
        """
        runs all validation techniques above
        :param params: list of entries for fem grid calculation
        :param values: list of results of params
        :param known_x: list of sampling points (matrix)
        :param known_fx: list of results for known_x
        :param vali_x: list of list of entry values
        :param vali_fx: list of results to each vali_x
        :param surro_func: pointer to the surrogates predict function
        :param surro_class: class of the used surrogate
        :param update_params: parameters to pass to surrogate on fitting
        :return:
        """
        res = ValidationResults()
        res.deviation = self.calc_deviation(params, values, surro_func)
        res.rmse = self.calc_rmse(vali_x, vali_fx, surro_func)
        res.mae = self.calc_mae(vali_x, vali_fx, surro_func)
        res.rae = self.calc_rae(vali_x, vali_fx, surro_func)
        res.press = self.calc_press(known_x, known_fx, surro_func, surro_class, update_params=update_params)
        return res


class ValidationResults():
    def __init__(self):
        self.deviation = 0.
        self.rmse = 0.
        self.mae = 0.
        self.rae = 0.
        self.press = 0.
