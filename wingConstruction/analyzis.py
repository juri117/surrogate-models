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
from scipy import interpolate

USED_CORES = Constants().config.getint('meta', 'used_cores')
NON_LINEAR = False

max_g = 2.5
safety_fac = 1.5
mtow = 27987.
fuel_mass_in_wings = 2*2659.
first_wing_struct_mass = 2*1000.
wing_load = ((mtow - fuel_mass_in_wings - first_wing_struct_mass) * 9.81) * max_g * safety_fac * 0.5
engine_weight = (873.1 + 251.9) * 9.81
engine_pos_y = 3.
wing_length = 12.87
chord_length = 3.
chord_height = 0.55

density = 2810 #kg/m^3
shear_strength = 5.72e8 #3.31e8 #Pa
max_shear_strength = shear_strength

element_size = 0.1


def new_project(project_name):
    #project_name = 'meshSize_r{:02d}_t{:5f}'.format(rib_count, shell_thick)
    pro = Project(project_name)
    pro.halfSpan = wing_length
    pro.boxDepth = chord_length * 0.4
    pro.boxHeight = chord_height
    pro.ribs = int(wing_length) + 1
    pro.enginePos = engine_pos_y
    pro.engineWeight = engine_weight
    pro.boxOverhang = 0.
    pro.forceTop = (2./3.) * wing_load
    pro.forceBot = (1./3.) * wing_load
    pro.elementSize = element_size
    # pro1.elementSize = 0.05
    pro.elemType = 'qu4'
    pro.shellThickness = 0.009
    return pro


def run_project(pro):
    pro.generate_geometry(nonlinear=NON_LINEAR)
    pro.solve()
    print('############ DONE ############')
    if not pro.errorFlag:
        if NON_LINEAR:
            pro.post_process(template='wing_post_nl_simple')
        else:
            pro.post_process(template='wing_post')
            #pro.postprocess(template='wing_post_max_mises_fixed')

    print('#########################################')
    print('finished: ' + pro.workingDir)
    return pro


def collect_results(pro):
    l = pro.validate_load('loadTop.frc')
    l += pro.validate_load('loadBot.frc')
    load_error = (-0.5*wing_load) - l
    if not pro.errorFlag:
        export_row = str(pro.elementSize) + ',' \
        + str(pro.geo.calc_span_division(pro.halfSpan)) + ',' \
        + str(pro.ribs) + ',' \
        + str(pro.shellThickness) + ',' \
        + str(pro.geo.calc_weight(density)) + ',' \
        + str(pro.clx.dispD3Min) + ','\
        + str(pro.clx.dispD3Max) + ','\
        + str(pro.clx.stressMisesMin) + ','\
        + str(pro.clx.stressMisesMax) + ',' \
        + str(load_error)+'\n'
        return export_row
    return ''


def pool_run(projects):
    start = time.time()
    with Pool(USED_CORES) as p:
        projects = p.map(run_project, projects)
    print("Time taken = {0:.5f}".format(time.time() - start))
    return projects


def main_run(cleanup=False):
    ribs = np.arange(6, 40, 2)
    ribs = list(ribs)
    thick = np.arange(0.006, 0.03, 0.002)
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
    output_f = open(Constants().WORKING_DIR + '/'
                   + output_file_name,
                   'w')
    output_f.write('elementSizes,spanElementCount,ribs,shellThickness,weight,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,loadError\n')

    for p in projects:
        outStr = collect_results(p)
        if outStr != '':
            output_f.write(outStr)
            output_f.flush()
        else:
            print('ERROR, empty data return')
    output_f.close()
    if cleanup:
        for p in projects:
            p.remove()
    print('DONE with ALL')
    return output_file_name


def read_data_file(output_file_name):
    file_path = Constants().WORKING_DIR + '/' + output_file_name
    data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
    ribs_raw = data[:, 2]
    shell_thick_raw = data[:, 3]
    weight_raw = data[:, 4]
    max_stress_raw = data[:, 8]
    max_disp_raw = data[:, 5]

    n_rib = len(set(ribs_raw))
    n_thick = len(set(shell_thick_raw))

    ribs = sorted(list(set(ribs_raw)))
    shell_thick = sorted(list(set(shell_thick_raw)))
    # shellThick = [x * 1000 for x in shellThick]
    weight = np.zeros((n_thick, n_rib))
    max_stress = np.zeros((n_thick, n_rib))
    max_disp = np.zeros((n_thick, n_rib))
    for i in range(0, len(ribs_raw)):
        weight[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = weight_raw[i]
        max_stress[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = max_stress_raw[i]
        max_disp[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = max_disp_raw[i]

    return ribs, shell_thick, max_stress, max_disp, weight

def plot_results(output_file_name):
    ribs, shell_thick, max_stress, max_disp, weight = read_data_file(output_file_name)
    n_rib = len(ribs)
    n_thick = len(shell_thick)
    opti_ribs = []
    opti_shell = []
    opti_weight = []
    rib_mat, shell_mat = np.meshgrid(ribs, shell_thick)
    # interpol weight
    f_weight = interpolate.interp2d(rib_mat, shell_mat, weight, kind='linear')

    plot1 = PlotHelper(['ribs', 'max stress'])
    for i in range(0, len(shell_thick)):
        stress = max_stress[i]
        #if np.min(stress) < max_shear_strength and np.max(stress) > max_shear_strength:
        #   optRibs = np.interp(max_shear_strength, stress, ribs)
        #   f = interp1d(stress, ribs, kind='linear')
        #   plot1.ax.plot([f(max_shear_strength)], [max_shear_strength], 'go')
        #   opti_ribs.append(f(max_shear_strength))
        #   opti_shell.append(shellThick[i])
        plot1.ax.plot(ribs, max_stress[i], label='shell= {:03f}'.format(shell_thick[i]))
    plot1.ax.plot(ribs, np.full((len(ribs), 1),max_shear_strength), 'r--', label='Limit-Load')
    plot1.finalize(legendNcol=2)
    #plot1.show()

    plot2 = PlotHelper(['shellthickness in mm', 'max stress'])
    for i in range(0, len(ribs)):
        stress = max_stress[:,i]
        plot2.ax.plot(shell_thick, stress, label='ribs= {:f}'.format(ribs[i]))
        if np.min(stress) < max_shear_strength and np.max(stress) > max_shear_strength:
            # optRibs = np.interp(max_shear_strength, stress, ribs)
            f = interp1d(stress, shell_thick, kind='linear')
            plot2.ax.plot([f(max_shear_strength)], [max_shear_strength], 'go')
            opti_ribs.append(ribs[i])
            opti_shell.append(f(max_shear_strength))
            opti_weight.append(f_weight(ribs[i], f(max_shear_strength))[0])
    plot2.ax.plot(shell_thick, np.full((len(shell_thick), 1), max_shear_strength), 'r--', label='Limit-Load')
    #plot2.ax.set_xlim((min([x * 1000 for x in shellThick]), max([x * 1000 for x in shellThick])))

    plot2.finalize(legendNcol=2)
    plot2.show()

    opti_min_index = opti_weight.index(min(opti_weight))
    print('opt wight {:03f}'.format(opti_weight[opti_min_index]))
    print('@ {:00f} ribs'.format(opti_ribs[opti_min_index]))
    print('@ {:05f} shell-thickness'.format(opti_shell[opti_min_index]))

    plot3d = PlotHelper(['ribs', 'shell thickness in m', 'mises stress'])
    max_stress[max_stress > 1.2*max_shear_strength] = np.nan
    color_map = plt.cm.jet(weight/np.max(weight))
    surf = plot3d.ax.plot_surface(rib_mat, shell_mat, max_stress, facecolors=color_map,
                    cstride=1,
                    rstride=1)
    optiLine = plot3d.ax.plot(opti_ribs, opti_shell, max_shear_strength, 'k--', label='Limit-Load')
    optiPoint = plot3d.ax.plot([opti_ribs[opti_min_index]], [opti_shell[opti_min_index]], max_shear_strength, 'ro', label='glob. optimum')
    m = cm.ScalarMappable(cmap=cm.jet)
    m.set_array(weight)
    cbar = plt.colorbar(m)
    cbar.set_label('structure weight in kg', rotation=270)
    limit = np.full((n_thick, n_rib),max_shear_strength)
    plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.1)
    plot3d.ax.set_zlim3d(0, max_shear_strength)
    plot3d.finalize()
    plot3d.show()


def convergence_analyzis_run(cleanup=False):
    sizes = np.arange(0.04, .26, 0.01)
    sizes = list(sizes)
    projects = []
    for s in sizes:
        project_name = 'meshSize_s{:5f}'.format(s)
        pro = new_project(project_name)
        pro.elementSize = s
        pro.ribs = 14
        projects.append(pro)

    projects = pool_run(projects)

    output_file_name = 'convAna_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
    outputF = open(Constants().WORKING_DIR + '/' + output_file_name, 'w')
    outputF.write(
        'elementSizes,spanElementCount,ribs,shellThickness,weight,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,loadError\n')

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
    #convergence_analyzis_run(cleanup=True)

    output_file_name = '2drun_2018-08-10_09_59_27.csv'
    output_file_name = main_run(cleanup=False)
    plot_results(output_file_name)