
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

# the smooth whole function
fx = np.linspace(0, 10, 1001)
fy = list(map(f,fx))

# now we pretend we only know a view points
px = [0., 2., 4., 6., 8., 10.]
py = list(map(f,px))
n = len(px)

#first fixed exponent here
p = 2.
#first fixed factor here
theta = 0.1
corMat = np.zeros((n,n))
diff = []
expFunc = []
for row in range(0, n):
    for column in range(0, n):
        corMat[row][column] = math.exp(-theta * abs(px[row] - px[column])**p)
        #since this is a symetric mat...
        #corMat[column][row] = corMat[row][column]
        diff.append(px[row] - px[column])
        expFunc.append(math.exp(-theta * abs(px[row] - px[column])**p))

fig, ax = plt.subplots()
#rc('text', usetex=True)
rc('font', **font)
rc('xtick', labelsize=FONT_SIZE)
rc('ytick', labelsize=FONT_SIZE)
plt.plot(diff, expFunc, 'bo', label='p=?, theta=?')
ax.legend(loc=1, ncol=1)#, mode="expand")
ax.set_xlabel('x_i - x', fontdict=font)
ax.set_ylabel('exp(-\theta |x_i-x|^p)', fontdict=font)
ax.tick_params(labelsize=16., length=6, width=2)
plt.tight_layout()
#plt.show()

standardDeviation = np.std(py)
sigma = standardDeviation
covMat = (sigma**2) * corMat.copy()
printMat(corMat, octave=True)

#U = np.linalg.cholesky(corMat)


y = np.array(py)
LnDetCorMat = np.log(np.linalg.det(corMat))
#LnDetPsi=2*sum(np.log(abs(np.diag(U))))
#mu=(np.transpose(np.ones((n,1)))*(U\(np.transpose(U)\y)))/(np.transpose(np.ones((n,1)))*(U\(np.transpose(U)\np.ones((n,1)))));
#SigmaSqr=(np.transpose((y-np.ones((n,1))*mu))*(U\(np.transpose(U)\(y-np.ones((n,1))*mu))))\n;
#NegLnLike=-1*(-(n/2)*np.log(SigmaSqr)-0.5*LnDetPsi);
one = np.ones((n,1)).flatten()
mu = (np.transpose(one)@np.linalg.inv(corMat)@y) / (np.transpose(one)@np.linalg.inv(corMat)@one)
SigmaSqr = (np.transpose(y - one*mu) @ np.linalg.inv(corMat) @ (y - one*mu)) / n

NegLnLike = (-1)*(-(n/2)*np.log(SigmaSqr) - 0.5*LnDetCorMat)
print('negLnLike = '+str(NegLnLike))




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


fig, ax = plt.subplots()
#rc('text', usetex=True)
rc('font', **font)
rc('xtick', labelsize=FONT_SIZE)
rc('ytick', labelsize=FONT_SIZE)
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
#plt.show()
