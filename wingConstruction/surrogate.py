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
from myLibs.Kriging import Kriging
from myLibs.Sampling import Sampling
from myLibs.Validation import Validation

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

# my guess
known_x_i = [0, 15,  0, 15, 0, 4, 5,  5, 8, 9, 11, 12, 14]
known_y_i = [0,  0, 13, 13, 5, 7, 3, 10, 2, 6, 11,  2,  5]
sample_indices = np.array([known_y_i, known_x_i]).T.tolist()


# latin hypercube
sam = Sampling()
sample_mat = sam.enhanced_latin_hypercube(14)
sample_indices = sam.bool_mat_to_list(sample_mat)


known_rib = []
known_shell = []
known_stress = []
for i in range(0, len(sample_indices)):
    known_rib.append(ribs[sample_indices[i][1]])
    known_shell.append(shell[sample_indices[i][0]])
    known_stress.append(stress[sample_indices[i][0]][sample_indices[i][1]])
known_params = [known_rib, known_shell]

print('sample plan using {:d} known values'.format(len(sample_indices)))

##################################################
# build surrogate model

krig = Kriging(known_params, known_stress)

##################################################
# fit surrogate model

#prev stored results:
krig.update_param([0.024898436078495855, 2482513.627234583], [1.8254148090594156, 2.0])

krig.optimize()

krig.plot_theta_likelihood_R2()
krig.plot_p_likelihood_R2()

minLike = krig.calc_likelihood()
print('minLike = ' + str(minLike))
print('@theta1 = ' + str(krig._theta[0]))
print('@theta2 = ' + str(krig._theta[1]))
print('@p1 = ' + str(krig._p[0]))
print('@p2 = ' + str(krig._p[1]))


##################################################
# validate

vali = Validation()
vali.calc_deviation(ribs, shell, stress, krig.predict)


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
realDat = plot3d.ax.plot_wireframe(rib_mat, shell_mat, stress, color='g', alpha=0.5, label='fem data')
#realDatMark = plot3d.ax.scatter(rib_mat, shell_mat, stress, c='g', marker='x', label='fem measurements')

#optiLine = plot3d.ax.plot(opti_ribs, opti_shell, max_shear_strength, 'k--', label='Limit-Load')
#optiPoint = plot3d.ax.plot([opti_ribs[opti_min_index]], [opti_shell[opti_min_index]], max_shear_strength, 'ro', label='glob. optimum')
#m = cm.ScalarMappable(cmap=cm.jet)
#m.set_array(weight)
#cbar = plt.colorbar(m)
#cbar.set_label('structure weight in kg', rotation=270)
limit = np.full((n_thick, n_rib),max_shear_strength)
plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.2, label='limit load')

ribs_sample = np.linspace(min(ribs), max(ribs), 200)
shell_sample = np.linspace(min(shell), max(shell), 200)
kritPlot = plot3d.plot_function_3D(krig.predict, ribs_sample, shell_sample, r'$f_{krig}$', color='b')
samplePoints = plot3d.ax.plot(known_rib, known_shell, known_stress, 'bo', label='sampling points')

plot3d.ax.set_zlim3d(0, max_shear_strength)
plot3d.finalize()
plot3d.show()

print('done')