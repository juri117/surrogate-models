__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :main file for testing
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import numpy as np
from datetime import datetime
from wingConstruction.utils.Constants import Constants
from wingConstruction.Project import Project
from utils.TimeTrack import TimeTrack
from utils.PlotHelper import PlotHelper

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm
from multiprocessing import Pool
import time
from scipy.interpolate import interp1d

USED_CORES = Constants().config.getint('meta', 'used_cores')
NON_LINEAR = False

mtow = 29257
fuel_mass_in_wings = 2*2659.
wing_load = (mtow - fuel_mass_in_wings) * 9.81
engine_mass = 1125 * 9.81
engine_pos_y = 3
wing_length = 12.87
chord_length = 3.
chord_height = 0.55

density = 2810 #kg/m^3
shear_strength = 3.31e8
max_g = 2.5
safety_fac = 1.5
max_shear_strength = shear_strength * max_g * safety_fac

element_size = 0.1


def new_project(project_name):
    #project_name = 'meshSize_r{:02d}_t{:5f}'.format(rib_count, shell_thick)
    pro = Project(project_name)
    pro.halfSpan = wing_length
    pro.boxDepth = chord_length * 0.4
    pro.boxHeight = chord_height
    pro.ribs = int(wing_length) + 1
    pro.boxOverhang = 0.
    pro.forceTop = -0.3 * wing_load
    pro.forceBot = -0.2 * wing_load
    pro.elementSize = element_size
    # pro1.elementSize = 0.05
    pro.elemType = 'qu8'
    pro.shellThickness = 0.009
    return pro


def run_project(pro):
    pro.generate_geometry(nonlinear=NON_LINEAR)
    pro.solve()
    print('############ DONE ############')
    if not pro.errorFlag:
        if NON_LINEAR:
            pro.postprocess(template='wing_post_nl_simple')
        else:
            pro.postprocess(template='wing_post_simple')
            pro.postprocess(template='wing_post_max_mises_fixed')

    print('#########################################')
    print('finished: ' + pro.workingDir)
    return pro


def collect_results(pro):
    l = pro.validate_load('loadTop.frc')
    l += pro.validate_load('loadBot.frc')
    loadError = (-0.5*wing_load) - l
    if not pro.errorFlag:
        exportRow = str(pro.elementSize) + ',' \
        + str(pro.ribs) + ',' \
        + str(pro.shellThickness) + ',' \
        + str(pro.geo.calc_weight(density)) + ',' \
        + str(pro.clx.dispD3Min) + ','\
        + str(pro.clx.dispD3Max) + ','\
        + str(pro.clx.stressMisesMin) + ','\
        + str(pro.clx.stressMisesMax) + ',' \
        + str(pro.clx.stressMisesMaxFixed) + ','\
        + str(pro.geo.calc_span_division(pro.halfSpan)) + ',' \
        + str(loadError)+'\n'
        return exportRow
    return ''

def pool_run(projects):
    start = time.time()
    with Pool(USED_CORES) as p:
        projects = p.map(run_project, projects)
    print("Time taken = {0:.5f}".format(time.time() - start))
    return projects

def main_run(cleanup=False):
    ribs = np.arange(1, 51, 2)
    ribs = list(ribs)
    thick = np.arange(0.0008, 0.0017, 0.00005)
    thick = list(thick)
    projects = []
    for r in ribs:
        for t in thick:
            project_name = 'pro_r{:02d}_t{:5f}'.format(r, t)
            pro = new_project(project_name)
            pro.ribs = r
            pro.shellThickness = t
            projects.append(pro)

    projects = pool_run(projects)

    output_file_name = '2drun_'+datetime.now().strftime('%Y-%m-%d_%H_%M_%S')+'.csv'
    outputF = open(Constants().WORKING_DIR + '/'
                   + output_file_name,
                   'w')
    outputF.write('elementSizes,ribs,shellThickness,weight,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,stressMisesMaxFixed,spanElementCount,loadError\n')

    for p in projects:
        outStr = collect_results(p)
        if outStr != '':
            outputF.write(outStr)
            outputF.flush()
        else:
            print('ERROR, empty data return')
    outputF.close()
    if cleanup:
        for p in projects:
            p.remove()
    print('DONE with ALL')
    return output_file_name


def plot_results(output_file_name):
    file_path = Constants().WORKING_DIR + '/' + output_file_name
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    ribsRaw = data[:, 1]
    shellThickRaw = data[:, 2]
    weightRaw = data[:, 3]
    maxStressRaw = data[:, 7]
    maxDispRaw = data[:, 4]

    nRib = len(set(ribsRaw))
    nThick = len(set(shellThickRaw))

    ribs = sorted(list(set(ribsRaw)))
    shellThick = sorted(list(set(shellThickRaw)))
    weight = np.zeros((nThick, nRib))
    maxStress = np.zeros((nThick, nRib))
    maxDisp = np.zeros((nThick, nRib))
    for i in range(0, len(ribsRaw)):
        weight[shellThick.index(shellThickRaw[i])][ribs.index(ribsRaw[i])] = weightRaw[i]
        maxStress[shellThick.index(shellThickRaw[i])][ribs.index(ribsRaw[i])] = maxStressRaw[i]
        maxDisp[shellThick.index(shellThickRaw[i])][ribs.index(ribsRaw[i])] = maxDispRaw[i]

    plot1 = PlotHelper(['ribs', 'max stress'])
    optiRibs = []
    optiShell = []
    for i in range(0, len(shellThick)):
        stress = maxStress[i]
        if np.min(stress) < max_shear_strength and np.max(stress) > max_shear_strength:
            #optRibs = np.interp(max_shear_strength, stress, ribs)
            f = interp1d(stress, ribs, kind='linear')
            plot1.ax.plot([f(max_shear_strength)], [max_shear_strength], 'go')
            optiRibs.append(f(max_shear_strength))
            optiShell.append(shellThick[i])

        plot1.ax.plot(ribs, maxStress[i], label='shell= {:03f}'.format(shellThick[i]))
    plot1.ax.plot(ribs, np.full((len(ribs), 1),max_shear_strength), 'r--', label='Limit-Load')
    plot1.finalize()
    #plot1.show()

    plot2 = PlotHelper(['ribs', 'max disp'])
    for i in range(0, len(shellThick)):
        plot2.ax.plot(ribs, maxStress[i], label='shell= {:03f}'.format(shellThick[i]))
    #plot1.ax.plot(ribs, np.full((len(ribs), 1), max_shear_strength), 'r--', label='Limit-Load')
    plot2.finalize()
    plot2.show()

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plotX, plotY = np.meshgrid(ribs, shellThick)
    maxStress[maxStress > 1.2*max_shear_strength] = np.nan
    #ax.plot_wireframe(plotX, plotY, maxStress, color='b', rstride=2, cstride=10)
    color_map = plt.cm.jet(weight/np.max(weight))
    surf = ax.plot_surface(plotX, plotY, maxStress, facecolors=color_map,
                    cstride=1,
                    rstride=1)
    optiLine = ax.plot(optiRibs, optiShell, max_shear_strength, 'k-')
    #cset = ax.contour(plotX, plotY, weight, 1000, zdir='z', offset=0, cmap=cm.coolwarm)
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(weight)
    plt.colorbar(m)
    limit = np.full((nThick, nRib),max_shear_strength)
    ax.plot_wireframe(plotX, plotY, limit, color='r', alpha=0.1)
    ax.set_zlim3d(0, max_shear_strength);
    plt.show()

def convergence_analyzis_run(cleanup=False):
    sizes = np.arange(0.04, .26, 0.01)
    sizes = list(sizes)
    projects = []
    for s in sizes:
        project_name = 'meshSize_s{:5f}'.format(s)
        pro = new_project(project_name)
        pro.elementSize = s
        projects.append(pro)

    projects = pool_run(projects)

    output_file_name = 'convAna_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    outputF = open(Constants().WORKING_DIR + '/' + output_file_name, 'w')
    outputF.write(
        'elementSizes,ribs,shellThickness,weight,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,stressMisesMaxFixed,spanElementCount,loadError\n')

    for p in projects:
        outStr = collect_results(p)
        if outStr != '':
            outputF.write(outStr)
            outputF.flush()
        else:
            print('ERROR, empty data return')
    outputF.close()
    if cleanup:
        for p in projects:
            p.remove()
    print('DONE with ALL')
    return output_file_name

if __name__ == '__main__':
    convergence_analyzis_run(cleanup=True)

    #output_file_name = 'convAna_2018-08-07_13_04_31.csv'
    #output_file_name = main_run(cleanup=True)
    #plot_results(output_file_name)