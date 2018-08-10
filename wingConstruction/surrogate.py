__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :creates the surrogate for wing-structure test
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
from datetime import datetime
from wingConstruction.utils.Constants import Constants
from wingConstruction.Project import Project
from wingConstruction.MultiRun import MultiRun
from utils.TimeTrack import TimeTrack
from utils.PlotHelper import PlotHelper

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from multiprocessing import Pool
import time
from scipy.interpolate import interp1d
from scipy import interpolate

RESULTS_FILE = '/dataOut/oldRun/2drun_2018-08-10_12_13_55.csv'

max_shear_strength = 5.72e8 * 2.5
multi = MultiRun()


##################################################
# collect data

ribs, shell, stress, disp, weight = multi.read_data_file(RESULTS_FILE)
n_rib = len(ribs)
n_thick = len(shell)
rib_mat, shell_mat = np.meshgrid(ribs, shell)

# input 0 is ribs
# input 1 is shell_thick

##################################################
# sample plan

known_x = [0, 15, 00, 15, 4, 5, 5, 8, 11, 14]
known_y = [0, 00, 13, 13, 7, 3, 10, 2, 11,  5]
known_rib = []
known_shell = []
known_stress = []
for i in range(0, len(known_x)):
    known_rib.append(ribs[known_x[i]])
    known_shell.append(shell[known_y[i]])
    known_stress.append(stress[known_y[i]][known_x[i]])


##################################################
# build surrogate model


##################################################
# optimize


##################################################
# plot it

plot3d = PlotHelper(['ribs', 'shell thickness in m', 'mises stress'])
#max_stress[max_stress > 1.2*max_shear_strength] = np.nan
#color_map = plt.cm.jet(weight/np.max(weight))
#surf = plot3d.ax.plot_surface(rib_mat, shell_mat, max_stress, facecolors=color_map,
#                cstride=1,
#                rstride=1)
realDat = plot3d.ax.plot_wireframe(rib_mat, shell_mat, stress, color='g', alpha=0.5, label='real data')

plot3d.ax.plot(known_rib, known_shell, known_stress, 'og', label='known points')

#optiLine = plot3d.ax.plot(opti_ribs, opti_shell, max_shear_strength, 'k--', label='Limit-Load')
#optiPoint = plot3d.ax.plot([opti_ribs[opti_min_index]], [opti_shell[opti_min_index]], max_shear_strength, 'ro', label='glob. optimum')
#m = cm.ScalarMappable(cmap=cm.jet)
#m.set_array(weight)
#cbar = plt.colorbar(m)
#cbar.set_label('structure weight in kg', rotation=270)
limit = np.full((n_thick, n_rib),max_shear_strength)
plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.1, label='limit load')
plot3d.ax.set_zlim3d(0, max_shear_strength)
plot3d.finalize()
plot3d.show()

print('done')