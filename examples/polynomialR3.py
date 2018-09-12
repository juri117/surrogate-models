
import numpy as np

from myLibs.Polynomial import Polynomial
from utils.PlotHelper import PlotHelper
from utils.TimeTrack import TimeTrack
from utils.samples import *

if __name__ == '__main__':
    t1 = TimeTrack('OverAllTimer')
    plt1 = PlotHelper(['Eingang 1', 'Eingang 2', 'Ausgang'], fancy=False, pgf=False)

    fx = np.linspace(-2, 12, 101)
    fy = np.linspace(-2, 12, 101)
    plt1.plot_function_3d(f_3D, fx, fy, r'$f_{original}$', color='r')
    # the smooth whole function

    # now we pretend we only know a view points
    pxEdge = [0., 2., 4., 6., 8., 10.]
    pyEdge = [0., 2., 4., 6., 8., 10.]
    px, py, pz = generate_sample_data(f_3D, pxEdge, pyEdge)
    knownParams = []
    knownParams.append(px)
    knownParams.append(py)

    scat1 = plt1.ax.scatter(px, py, pz, c='r', marker='o', s=10, label=r'St\"utzstellen')

    poly = Polynomial(knownParams, pz)
    color = ['dodgerblue', 'gold', 'fuchsia', 'lime', 'lightsalmon', 'maroon']
    ic = 0
    for o in [6]:
        poly.update_param(o)
        plt1.plot_function_3d(poly.predict,
                              fx,
                              fy,
                              r'$\widehat{f}_{Poly} mit o = ' + str(o) + '$',
                              color=color[ic])
        poly.generate_formula()
        ic += 1

    t1.toc()

    plt1.ax.view_init(20, 50)
    plt1.ax.set_zlim3d(np.min(pz), max(pz) * 1.2)
    plt1.finalize(width=8, height=5)
    plt1.show()


