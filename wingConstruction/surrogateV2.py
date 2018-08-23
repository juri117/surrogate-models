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
from wingConstruction.fem.WingConstructionV4 import WingConstruction

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from multiprocessing import Pool
import time
from scipy.interpolate import interp1d
from scipy import interpolate
from scipy.optimize import minimize
from scipy import optimize

RESULTS_FILE = '/2drun_2018-08-23_16_49_18_final01_cruiseLoad.csv'

wing_length = 12.87
chord_length = 3.
chord_height = 0.55
density = 2810 #kg/m^3
max_shear_strength = 5.72e8 / 1.5
multi = MultiRun()


##################################################
# collect data

ribs, shell, stress, disp, weight = multi.read_data_file(RESULTS_FILE)
n_rib = len(ribs)
n_thick = len(shell)
rib_mat, shell_mat = np.meshgrid(ribs, shell)

# input 0 is ribs
# input 1 is shell_thick

# reduce data range
i_rib_0 = 0
i_rib_n = 15
i_shell_0 = 0
i_shell_n = 15
newRib = []
newShell = []
newStress = np.zeros((i_shell_n-i_shell_0, i_rib_n-i_rib_0))
newWeight = np.zeros((i_shell_n-i_shell_0, i_rib_n-i_rib_0))
for r in range(0, i_rib_n-i_rib_0):
    newRib.append(ribs[i_rib_0+r])
    for s in range(0, i_shell_n-i_shell_0):
        newStress[s][r] = stress[i_shell_0+s][i_rib_0+r]
        newWeight[s][r] = weight[i_shell_0 + s][i_rib_0 + r]

for s in range(0, i_shell_n-i_shell_0):
    newShell.append(shell[i_shell_0+s])

ribs = newRib
shell = newShell
stress = newStress
weight = newWeight

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
krig.update_param([0.002261264770141511, 277826.21903867245], [1.8766170168043503, 1.9959876593551822])

#krig.optimize()

#krig.plot_theta_likelihood_R2()
#krig.plot_p_likelihood_R2()

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

def shell_predict(shell_thick, krig_inst, rib_num):
    #krig_inst = args[0]
    #rib_num = args[1]
    stress_val = krig_inst.predict([rib_num, shell_thick])
    #return stress_val
    return stress_val - max_shear_strength

opti_ribs = []
opti_shell = []
opti_stress = []
opti_weights = []
used_ribs_indices = list(range(ribs.index(min(known_rib)), ribs.index(max(known_rib))))
used_shell_indices = list(range(shell.index(min(known_shell)), shell.index(max(known_shell))))
for i in used_ribs_indices:
    #stress = stress[:,i]
    # SLSQP: proplem; find local min not glob. depending on init-vals
    init_guess = shell[int(len(known_shell)/2.)]
    bnds = [(min(known_shell), max(known_shell))]
    #res = minimize(shell_predict, init_guess, args=[krig, ribs[i]], method='SLSQP', tol=1e-6, options={'disp': True, 'maxiter': 99999}, bounds=bnds)
    #opti_shell.append(res.x[0])
    root = optimize.newton(shell_predict, init_guess, args=[krig, ribs[i]])
    opti_ribs.append(ribs[i])
    opti_shell.append(root)
    opti_stress.append(krig.predict([ribs[i], root]))
    weight = WingConstruction.calc_weight_stat(wing_length,
                                               chord_length,
                                               chord_height,
                                               ribs[i],
                                               root,
                                               density)
    opti_weights.append(weight)
    """
    plt.clf()
    sll = np.linspace(min(shell), max(shell), 500)
    strs = []
    for s in sll:
        strs.append(shell_predict(s, krig, ribs[i]))
    plt.plot(sll, strs, 'b-')
    opt_strs = shell_predict(root, krig, ribs[i])
    plt.plot([root], [opt_strs], 'ro')
    plt.grid(True)
    plt.show()
    """

# exclude model edges from opti vals
opti_ribs = opti_ribs[:-1]
opti_ribs = opti_ribs[1:]
opti_shell = opti_shell[:-1]
opti_shell = opti_shell[1:]
opti_stress = opti_stress[:-1]
opti_stress = opti_stress[1:]
opti_weights = opti_weights[:-1]
opti_weights = opti_weights[1:]

best_i = opti_weights.index(min(opti_weights))

optWeightPlot = PlotHelper(['ribs', 'weight'])
optWeightPlot.ax.plot(opti_ribs, opti_weights, 'b-')
optWeightPlot.ax.plot([opti_ribs[best_i]], opti_weights[best_i], 'rx', label='minimum')
optWeightPlot.finalize()

##################################################
# plot it

# convert all to np arrays
shell = np.array(shell)
opti_shell = np.array(opti_shell)
known_shell = np.array(known_shell)

plot3d = PlotHelper(['ribs', 'shell thickness in mm', 'mises stress'], fancy=False, font_size=16)

#realDat = plot3d.ax.plot_wireframe(rib_mat, shell_mat, stress, color='g', alpha=0.5, label='fem data')

# plot FEM data as lines
# plot FEM data as lines
for i in range(0,len(ribs)):
    if i == 0:
        plot3d.ax.plot(np.ones((len(shell)))*ribs[i], shell*1000., stress[:,i], 'g-', lw=3., label='fem data')
    else:
        plot3d.ax.plot(np.ones((len(shell))) * ribs[i], shell*1000., stress[:, i], 'g-', lw=3.)

#realDatMark = plot3d.ax.scatter(rib_mat, shell_mat, stress, c='g', marker='x', label='fem measurements')

# plot limit load as wireframe
#limit = np.full((n_thick, n_rib),max_shear_strength)
#plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.2, label='limit load')

# plot surrogate model as wireframe
ribs_sample = np.linspace(min(known_rib), max(known_rib), 200)
shell_sample = np.linspace(min(known_shell), max(known_shell), 200)
kritPlot = plot3d.plot_function_3D(krig.predict, ribs_sample, shell_sample, r'$f_{krig}$', color='b', scale=[1.,1000.,1.])
samplePoints = plot3d.ax.plot(known_rib, known_shell*1000., known_stress, 'bo', label='sampling points')

# plot limit load line
plot3d.ax.plot(opti_ribs, opti_shell*1000., opti_stress, 'k--', lw=3., label='max. stress line')
# plot optimal point
plot3d.ax.plot([opti_ribs[best_i]], [opti_shell[best_i]*1000.], [opti_stress[best_i]], 'rx', markersize=12, markeredgewidth=5, label='global optimum')
plot3d.ax.locator_params(nbins=7, axis='y')

plot3d.ax.set_zlim3d(np.min(np.array(stress)), max_shear_strength*1.2)

plot3d.finalize(height=7, width=9, legendLoc=8, legendNcol=3, bbox_to_anchor=(0.5, -0.0), tighten_layout=True)
plot3d.ax.view_init(18, 40)
plot3d.save('../dataOut/surroPlotV2.pdf')
plot3d.show()

print('done')