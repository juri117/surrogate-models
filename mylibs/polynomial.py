__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Polynomial Model
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import numpy as np

VARS = ['x', 'y', 'z']


class Polynomial:

    def __init__(self, known_in, known_val):
        """
        :param known_in: list of lists with input sample points
        :param known_val: list of results for the known_in
        """
        self._known_in = np.array(known_in)
        self._known_val = np.array(known_val)
        if len(self._known_in.shape) == 1:
            self._known_in = self._known_in.reshape((self._known_in.shape[0], 1))
        self._k = self._known_in.shape[1]
        self._n = self._known_in.shape[0]
        self._order = 2

    def train(self):
        """
        trains the surrogate if available
        :return: None
        """
        self.update_param(self._order)

    def update_param(self, order):
        """
        updates the parameters of the surrogate model
        :param order: the maximum polynomial order
        :return: None
        """
        self._order = order
        self._calc_vandermonde_mat()
        self._calc_weights()

    def _calc_term_count(self):
        iw = 1
        for o in range(1, self._order + 1):
            for ik in range(0, self._k):
                iw += 1
                for ikc in range(0, self._k):
                    if ikc > ik and 2 * o < self._order + 1:
                        iw += 1
                    if ikc != ik:
                        for ioc in range(1, min(o, (self._order + 1) - o)):
                            iw += 1
        return iw

    def _calc_vandermonde_mat(self):
        vander = np.zeros((self._n, self._calc_term_count()))
        for i in range(0, self._n):
            iw = 0
            vander[i][iw] = 1
            iw += 1
            for o in range(1, self._order+1):
                for ik in range(0, self._k):
                    vander[i][iw] = self._known_in[i][ik] ** o
                    iw += 1
                    for ikc in range(0, self._k):
                        if ikc > ik and 2*o < self._order+1:
                            vander[i][iw] = (self._known_in[i][ik] ** o) * (self._known_in[i][ikc] ** o)
                            iw += 1
                        if ikc != ik:
                            for ioc in range(1, min(o, (self._order+1)-o)):
                                vander[i][iw] = (self._known_in[i][ik] ** o) * (self._known_in[i][ikc] ** ioc)
                                iw += 1
        # delete unused columns
        for i in range(0, self._calc_term_count() - iw):
            vander = np.delete(vander, -1, axis=1)
            print('KICK OUT')
        self._vander = vander

    def _calc_weights(self):
        # moore-penrose pseudo-inverse
        pin_vander = np.linalg.pinv(self._vander)
        weights = pin_vander @ self._known_val
        self._weights = weights

    def predict(self, x_pred):
        """
        predicts a value from the surrogate model
        :param x_pred: vector of input values
        :return: result value
        """
        fx = 0.
        iw = 0
        fx += self._weights[iw]
        iw += 1
        for o in range(1, self._order + 1):
            for ik in range(0, self._k):
                fx += self._weights[iw] * x_pred[ik] ** o
                iw += 1
                for ikc in range(0, self._k):
                    if ikc > ik and 2 * o < self._order + 1:
                        fx += self._weights[iw] * (x_pred[ik]**o) * (x_pred[ikc]**o)
                        iw += 1
                    if ikc != ik:
                        for ioc in range(1, min(o, (self._order + 1) - o)):
                            fx += self._weights[iw] * (x_pred[ik]**o) * (x_pred[ikc]**ioc)
                            iw += 1
        return fx

    def generate_formula(self):
        str_print = ''
        iw = 0
        str_print += '{:f}'.format(self._weights[iw])
        iw += 1
        for o in range(1, self._order + 1):
            for ik in range(0, self._k):
                str_print += ' + {:f} * {:s}^({:d})'.format(self._weights[iw], VARS[ik], o)
                iw += 1
                for ikc in range(0, self._k):
                    if ikc > ik and 2 * o < self._order + 1:
                        str_print += ' + {:f} * {:s}^({:d}) * {:s}^({:d})'.format(self._weights[iw], VARS[ik], o, VARS[ikc], o)
                        iw += 1
                    if ikc != ik:
                        for ioc in range(1, min(o, (self._order + 1) - o)):
                            str_print += ' + {:f} * {:s}^({:d}) * {:s}^({:d})'.format(self._weights[iw], VARS[ik], o, VARS[ikc],
                                                                            ioc)
                            iw += 1
        str_print = str_print.replace('+ -', '- ')
        print(str_print)
        return str_print

    def get_order(self):
        return self._order

    def get_weights(self):
        return self._weights
