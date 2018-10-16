
import numpy as np

from mylibs.rbf import RBF
from myutils.plot_helper import PlotHelper
from myutils.time_track import TimeTrack
from myutils.samples import *

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
    knownParams = np.array([px, py]).T

    scat1 = plt1.ax.scatter(px, py, pz, c='r', marker='o', s=10, label=r'St\"utzstellen')

    rbf = RBF(knownParams, pz)
    a = 0.17
    rbf.update_param(a, 'gaus')

    plt1.plot_function_3d(rbf.predict, fx, fy, r'$\widehat{f}_{RBF}$', color='b')
    t1.toc()

    plt1.ax.view_init(20, 50)
    plt1.finalize(width=8, height=5)
    plt1.show()


