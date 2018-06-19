
import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib import rc
import matplotlib
import scipy

FONT_SIZE=14
font = {'family':'sans-serif', 'size':FONT_SIZE}

def printMat(mat, octave=False):
    if octave:
        print('['+';\n'.join([''.join(['{:10.3} '.format(item) for item in row]) for row in mat])+']')
    else:
        print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4

#calcs the correlation matrix
def calc_corMat(x, theta, p, plot_it=False):
    n = len(x)
    corMat = np.zeros((n, n))
    diff = []
    expFunc = []
    for row in range(0, n):
        for column in range(0, n):
            corMat[row][column] = math.exp(-theta * abs(x[row] - x[column]) ** p)
            # since this is a symetric mat...
            # corMat[column][row] = corMat[row][column]
            diff.append(x[row] - x[column])
            expFunc.append(math.exp(-theta * abs(x[row] - x[column]) ** p))

    if plot_it:
        fig, ax = plt.subplots()
        # rc('text', usetex=True)
        rc('font', **font)
        rc('xtick', labelsize=FONT_SIZE)
        rc('ytick', labelsize=FONT_SIZE)
        plt.plot(diff, expFunc, 'bo', label='p=?, theta=?')
        ax.legend(loc=1, ncol=1)  # , mode="expand")
        ax.set_xlabel('x_i - x', fontdict=font)
        ax.set_ylabel('exp(-\theta |x_i-x|^p)', fontdict=font)
        ax.tick_params(labelsize=16., length=6, width=2)
        plt.tight_layout()
        plt.show()
    return corMat

#calcs the likelyhood
def calc_likelihood(x, y, theta, p):
    corMat = calc_corMat(x, theta, p)
    n = len(y)
    y = np.array(y)
    LnDetCorMat = np.log(np.linalg.det(corMat))
    one = np.ones((n, 1)).flatten()
    mu = (np.transpose(one) @ np.linalg.inv(corMat) @ y) / (np.transpose(one) @ np.linalg.inv(corMat) @ one)
    SigmaSqr = (np.transpose(y - one * mu) @ np.linalg.inv(corMat) @ (y - one * mu)) / n
    NegLnLike = (-1) * (-(n / 2) * np.log(SigmaSqr) - 0.5 * LnDetCorMat)
    return NegLnLike

def predict(x_pred, x, y, theta, p):
    n = len(y)
    corMat = calc_corMat(x, theta, p)
    one = np.ones((n, 1)).flatten()
    mu = (np.transpose(one) @ np.linalg.inv(corMat) @ y) / (np.transpose(one) @ np.linalg.inv(corMat) @ one)
    psi = np.ones((n, 1)).flatten()
    for i in range(0, len(psi)):
        psi[i] = math.exp(-theta * abs(x[i] - x_pred) ** p)
    fx = mu + np.transpose(psi) @ np.linalg.inv(corMat) @ (y - one * mu)
    return fx

if __name__ == '__main__':
    # the smooth whole function
    fx = np.linspace(0, 10, 1001)
    fy = list(map(f,fx))

    # now we pretend we only know a view points
    px = [0., 2., 4., 6., 8., 10.]
    py = list(map(f,px))


    #first fixed exponent here
    p = 2.
    #first fixed factor here
    theta = 1.


    standardDeviation = np.std(py)
    sigma = standardDeviation
    #covMat = (sigma**2) * corMat.copy()
    #printMat(corMat, octave=True)

    #U = np.linalg.cholesky(corMat)

    NegLnLike = calc_likelihood(px, py, theta, p)
    print('negLnLike = ' + str(NegLnLike))

    thetas = np.linspace(0.01, 1000, 1000)
    likely = []
    for thet in thetas:
        likely.append(calc_likelihood(px, py, thet, 2))
    #plt.plot(thetas, likely)
    #plt.show()

    fig, ax = plt.subplots()
    # rc('text', usetex=True)
    rc('font', **font)
    rc('xtick', labelsize=FONT_SIZE)
    rc('ytick', labelsize=FONT_SIZE)


    #fx = np.linspace(0., 10., 100+1)
    krigY = np.zeros((len(fx),1)).flatten()
    for i in range(0, len(fx)):
        krigY[i] = predict(fx[i], px, py, theta, p)
    ax.plot(fx, krigY, 'b-', label='prediction')
    #plt.show()





    """
    dist = []
    diff = []
    for i in range(0, len(px)):
        for i2 in range(i, len(px)):
            dist.append(abs(px[i] - px[i2]))
            diff.append(0.5 * (py[i] - py[i2])**2)
    
    plt.plot(dist, diff, 'bo')
    plt.show()
    """



    ax.plot(fx, fy, 'r-', label=r'$f_{original}$')
    ax.plot(px, py, 'ro', label=r'St\"utzstellen', markersize=10)
    ax.legend(loc=3, ncol=2, mode="expand")
    ax.set_xlabel('Eingang', fontdict=font)
    ax.set_ylabel('Ausgang', fontdict=font)
    ax.tick_params(labelsize=16., length=6, width=2)
    fig.set_size_inches(8, 5)
    plt.tight_layout()
    plt.savefig('dataOut/radialBasisR2.svg')
    plt.savefig('dataOut/radialBasisR2.pdf')
    plt.show()
