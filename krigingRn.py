
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
import matplotlib
import scipy
from scipy.optimize import minimize

FONT_SIZE=14
font = {'family':'sans-serif', 'size':FONT_SIZE}

def printMat(mat, octave=False):
    if octave:
        print('['+';\n'.join([''.join(['{:10.3} '.format(item) for item in row]) for row in mat])+']')
    else:
        print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f(x, y):
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4 + 0.05*y**2 - 0.001*y**4 - 0.005*y**3


#calcs the correlation matrix
def calc_corMat(known_x, theta, p, plot_it=False):
    k = known_x.shape[0]
    n = known_x.shape[1]
    corMat = np.zeros((n, n))
    diff = []
    expFunc = []
    for row in range(0, n):
        for column in range(0, n):
            sum = 0.
            for ik in range(0, k):
                sum += theta[ik] * (abs(known_x[ik][row] - known_x[ik][column])**p[ik])
            corMat[row][column] = math.exp(-sum)
            # since this is a symetric mat...
            #corMat[column][row] = corMat[row][column]
            #diff.append(knwon_x[row] - knwon_x[column])
            #expFunc.append(math.exp(-sum))
    """
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
    """
    return corMat

#calcs the likelyhood
def calc_likelihood(known_x, known_val, theta, p):
    corMat = calc_corMat(known_x, theta, p)
    k = known_x.shape[0]
    n = known_x.shape[1]
    known_val = np.array(known_val)
    #LnDetCorMat = np.log(np.linalg.det(corMat))
    LnDetCorMat = np.linalg.slogdet(corMat)[1]
    one = np.ones((n, 1)).flatten()
    mu = (np.transpose(one) @ np.linalg.inv(corMat) @ known_val) / (np.transpose(one) @ np.linalg.inv(corMat) @ one)
    SigmaSqr = (np.transpose(known_val - one * mu) @ np.linalg.inv(corMat) @ (known_val - one * mu)) / n
    NegLnLike = (-1) * (-(n / 2) * np.log(SigmaSqr) - 0.5 * LnDetCorMat)
    if NegLnLike == float('nan'):
        print('Error: nan')
    return NegLnLike

def calc_likelihood_opti(params, *args):
    NegLnLike = calc_likelihood(args[0], args[1], params, args[2])
    print(str(NegLnLike))
    return NegLnLike

def predict(x_pred, known_x, knwon_val, theta, p):
    k = known_x.shape[0]
    n = known_x.shape[1]
    corMat = calc_corMat(known_x, theta, p)
    one = np.ones((n, 1)).flatten()
    mu = (np.transpose(one) @ np.linalg.inv(corMat) @ knwon_val) / (np.transpose(one) @ np.linalg.inv(corMat) @ one)
    psi = np.ones((n, 1)).flatten()
    for i in range(0, n):
        sum = 0.
        for ik in range(0, k):
            sum += theta[ik] * (abs(known_x[ik][i] - x_pred[ik]) ** p[ik])
            psi[i] = math.exp(-sum)

        #psi[i] = math.exp(-theta * abs(knwon_x[i] - x_pred) ** p)
    fx = mu + np.transpose(psi) @ np.linalg.inv(corMat) @ (knwon_val - one * mu)
    return fx

if __name__ == '__main__':
    # the smooth whole function
    fx = np.linspace(-2, 12, 101)
    fy = np.linspace(-2, 12, 101)
    fz = np.zeros((len(fx), len(fy)))
    for iX in range(0,len(fx)):
        for iY in range(0,len(fy)):
            fz[iX][iY] = f(fx[iX], fy[iY])

    # now we pretend we only know a view points
    knownParams = []
    pxEdge = [0., 2., 4., 6., 8., 10.]
    pyEdge = [0., 2., 4., 6., 8., 10.]
    px = []
    py = []
    pz = []
    for iX in range(0, len(pxEdge)):
        for iY in range(0, len(pyEdge)):
            px.append(pxEdge[iX])
            py.append(pyEdge[iY])
            pz.append(f(pxEdge[iX], pyEdge[iY]))
    knownParams.append(px)
    knownParams.append(py)
    knownParams = np.array(knownParams)
    knownValues = np.array(pz)

    p = [2., 2.]

    #thetas = np.linspace(0.001, 0.1, 100 + 1)
    thetas = np.logspace(-3, 1, num=20)
    likelyX1 = []
    likelyX2 = []
    likely = np.zeros((len(thetas), len(thetas)))
    for i1 in range(0, len(thetas)):
        print(str(i1) + ' / ' + str(len(thetas)))
        for i2 in range(0, len(thetas)):
            #likelyX1.append(calc_likelihood(knownParams, knownValues, [thetas[i1], thetas[i2]], p))
            #likelyX2.append(calc_likelihood(knownParams, knownValues, [1., thet], p))
            likely[i1][i2] = calc_likelihood(knownParams, knownValues, [thetas[i1], thetas[i2]], p)

    x0 = [1., 1.]
    bnds = [(0.0001, 1000.), (0.0001, 1000.)]
    opt = {}
    opt['disp'] = True
    opt['maxiter'] = 99999
    res = minimize(calc_likelihood_opti, x0, args=(knownParams, knownValues, p), method='SLSQP', tol=1e-6, options=opt, bounds=bnds)

    bestThetaX1 = res.x[0]
    minLikeX1 = calc_likelihood(knownParams, knownValues, [bestThetaX1, 1.], p)
    bestThetaX2 = res.x[1]
    minLikeX2 = calc_likelihood(knownParams, knownValues, [1., bestThetaX2], p)
    bestTheta = [res.x[0], res.x[1]]
    minLike = calc_likelihood(knownParams, knownValues, [bestThetaX1, bestThetaX2], p)
    print('minLike = ' + str(minLike))
    print('@theta1 = ' + str(bestThetaX1))
    print('@theta2 = ' + str(bestThetaX2))
    fig, ax = plt.subplots()
    ax.set_yscale('log')
    ax.set_xscale('log')
    pcol = ax.pcolor(thetas, thetas, likely)
    fig.colorbar(pcol)
    ax.plot(bestThetaX1, bestThetaX2, 'rx')
    #ax.plot([bestThetaX1], [bestThetaX2], 'rx')

    #plt.semilogx(thetas, likelyX1, '-b')
    #plt.semilogx(thetas, likelyX2, '-g')
    #plt.semilogx(bestThetaX1, minLikeX1, 'bx')
    #plt.semilogx(bestThetaX2, minLikeX2, 'gx')
    #plt.semilogx([0, 1], [minLike, minLike], 'rx')
    plt.show()



    #plot original function and points
    plotX, plotY = np.meshgrid(fx, fy)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    # rc('text', usetex=True)
    rc('font', **font)
    surf1 = ax.plot_wireframe(plotX, plotY, fz, color='r', label=r'$f_{original}$', rcount=20, ccount=20, linewidths=1,
                              alpha=0.5)  # , rstride=1, cstride=1)#, cmap=cm.coolwarm) # ,linewidth=0, antialiased=False

    scat1 = ax.scatter(py, px, pz, c='r', marker='o', s=10, label=r'St\"utzstellen')





    krigingSol = np.zeros((len(fx), len(fy)))
    for iX in range(0, len(fx)):
        print(str(iX) + ' of ' + str(len(fx)))
        for iY in range(0, len(fy)):
            coords = [fx[iX], fy[iY]]
            krigingSol[iX][iY] = predict(coords, knownParams, knownValues, bestTheta, p)

    surf2 = ax.plot_wireframe(plotX, plotY, krigingSol, color='b', label=r'$f_{kriging}$', rcount=20, ccount=20, linewidths=1,
                              alpha=0.5)  # , rstride=1, cstride=1)#, cmap=cm.coolwarm) # ,linewidth=0, antialiased=False

    ax.view_init(20, 50)
    rc('xtick', labelsize=16)
    rc('ytick', labelsize=16)
    ax.set_xlabel('Eingang 1', fontdict=font)
    ax.set_ylabel('Eingang 2', fontdict=font)
    ax.set_zlabel('Ausgang', fontdict=font)
    ax.xaxis._axinfo['label']['space_factor'] = 4
    ax.tick_params(labelsize=16., length=6, width=2)
    fig.set_size_inches(8, 5)
    # plt.tight_layout()
    ax.legend()
    ax.autoscale_view(tight=True)
    plt.savefig('dataOut/krigingRn.svg')
    plt.savefig('dataOut/krigingRn.pdf')

    # for angle in range(0, 360):
    #    ax.view_init(30, angle)
    #    plt.draw()
    #    print(str(angle))
    #    plt.pause(.001)
    plt.show()
    print('done')














