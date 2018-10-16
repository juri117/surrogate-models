import math

def printMat(mat):
    print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))

# this is our sample function
def f_3D(coords):
    x = coords[0]
    y = coords[1]
    return math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4 + 0.05*y**2 - 0.001*y**4 - 0.005*y**3

def f_2D(x):
    return 6 - (math.sin(x) + 0.95 + 0.075*x**2 - 0.001*x**4)

def f_2D_big(x):
    return math.sin(x)+3.+0.035*x**2+0.003*x**3-0.00045*x**4+0.0000005*x**6


def generate_sample_data(f_3D, pxEdge, pyEdge):
    px = []
    py = []
    pz = []
    for iX in range(0, len(pxEdge)):
        for iY in range(0, len(pyEdge)):
            px.append(pxEdge[iX])
            py.append(pyEdge[iY])
            pz.append(f_3D([pxEdge[iX], pyEdge[iY]]))
    return px, py, pz

def printMat(mat, octave=False):
    if octave:
        print('['+';\n'.join([''.join(['{:10.3} '.format(item) for item in row]) for row in mat])+']')
    else:
        print('\n'.join([''.join(['{:10.3}\t'.format(item) for item in row]) for row in mat]))