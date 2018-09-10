__author__ = "Juri Bieler"
__version__ = "0.0.1"
__email__ = "juribieler@gmail.com"
__status__ = "Development"

# ==============================================================================
# description     :n-dimensional Sampling plans
# date            :2018-07-23
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================


import numpy as np
import math
import sys


class Hammersley:

    def __init__(self):
        print('done')

    '''
    :param n
    :return the n-th prime number (starting from n=0 => 2)
    '''
    def prime(self, n):
        n += 1
        prime_list = [2]
        num = 3
        while len(prime_list) < n:
            for p in prime_list:
                if num % p == 0:
                    break
            else:
                prime_list.append(num)
            num += 2
        return prime_list[-1]

    def hammersley(self, i, m, n):
        #  Parameters:
        #    Input, integer I, the index of the element of the sequence.
        #    0 <= I.
        #    Input, integer M, the spatial dimension.
        #    1 <= M <= 100.
        #    Input, integer N, the "base" for the first component.
        #    1 <= N.
        #    Output, real R(M), the element of the sequence with index I.
        import numpy as np
        i = int(i)
        t = np.ones(m - 1)
        t = i * t
        # Carry out the computation.
        prime_inv = np.zeros(m - 1)
        for j in range(0, m - 1):
            prime_inv[j] = 1.0 / float(self.prime(j))

        r = np.zeros(m)
        r[0] = float(i % (n + 1)) / float(n)

        while (0 < np.sum(t)):
            for j in range(0, m - 1):
                d = (t[j] % self.prime(j))
                r[j + 1] = r[j + 1] + float(d) * prime_inv[j]
                prime_inv[j] = prime_inv[j] / self.prime(j)
                t[j] = (t[j] // self.prime(j))
        return r

    def generate_sample_plan(self, point_count, dimension, bounds, base=None):
        usedBase = point_count
        if not base is None:
            usedBase = base
        points = []
        for i in range(0, point_count):
            normPoint = self.hammersley(i, dimension, usedBase)
            scaledPoint = []
            for d in range(0, dimension):
                scaledPoint.append(bounds[d][0] + (normPoint[d] * (bounds[d][1]- bounds[d][0])))
            points.append(scaledPoint)
        return points



if __name__ == '__main__':
    import matplotlib.pyplot as plt

    ham = Hammersley()

    #sys.exit(15)


    for i in range(0, 100):
        point = ham.hammersley(i, 2, 100)
        plt.plot([point[0]], [point[1]], 'bo')
    plt.show()

    samples = ham.generate_sample_plan(14, 2, [(5, 20), (0.01, 0.05)])
    for i in range(0, 14):
        plt.plot([samples[i][0]], [samples[i][1]], 'bo')
    plt.show()