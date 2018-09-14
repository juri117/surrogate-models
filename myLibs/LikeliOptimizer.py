__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Optimizer if the global min is hard to find
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

from scipy.optimize import minimize
import numpy as np

class LikeliOptimizer:

    def __init__(self, debug=False):
        self.debug = debug
        self.pBounds = (1., 2.)
        self.thetaBoundsExp = (-5, 10)

    def find(self, func, dimensions):
        opt={'disp': False, 'maxiter': 5e6}
        guess = self.generate_grid(func, dimensions, 15, 5)

        bnds = []
        for i in range(0, dimensions):
            bnds.append((-5, +10))
            #bnds.append((max(1e-5, guess[i] * 0.1), min(1e+10, guess[i] * 10)))
        for i in range(0, dimensions):
            bnds.append((1., 2.))
            #bnds.append((max(1., guess[dimensions+i] * 0.5), min(2., guess[dimensions+i] * 1.5)))

        minima_res = minimize(func, guess, method='SLSQP', tol=1e-8,
                       options=opt, bounds=bnds)
        return minima_res

    def generate_grid(self, func, dimensions, theta_sections, p_sections):
        minima_val = float('inf')
        minima_param = None
        ps = np.linspace(1., 2., num=p_sections)
        #thetas = np.logspace(-5, 9, num=theta_sections)
        thetas = np.linspace(-5, 9, num=theta_sections)
        p_param_i = np.zeros((dimensions))
        t_param_i = np.zeros((dimensions))

        params = []

        for ip in range(len(ps)**dimensions):
            #if p_param_i[ip] >= len(ps)-1:
            #    p_param_i[ip] = 0
            #else:
            #    p_param_i[ip] += 1
            t_param_i = np.zeros((dimensions))
            for it in range(len(thetas)**dimensions):
                #if t_param_i[it] >= len(thetas) - 1:
                #    t_param_i[it] = 0
                #else:
                #    t_param_i[it] += 1
                # run calc
                param = []
                for i in range(0, dimensions):
                    param.append(thetas[int(t_param_i[i])])
                for i in range(0, dimensions):
                    param.append(ps[int(p_param_i[i])])
                #print(str(param))
                res_val = func(param)
                # compare if new min
                if res_val < minima_val:
                    minima_val = res_val
                    minima_param = param
                # increment theta
                t_param_i = self.increase_i(t_param_i, len(thetas))
                params.append(param)
            # increment p
            p_param_i = self.increase_i(p_param_i, len(ps))
        return minima_param


    def increase_i(self, i_list, max_len):
        for i in range(0, len(i_list)):
            # -1 here because we start indexing at 0
            if i_list[i] >= max_len - 1:
                i_list[i] = 0
            else:
                i_list[i] += 1
                return i_list
        return i_list
