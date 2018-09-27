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
from wingConstruction.wingUtils.Constants import Constants
from wingConstruction.Project import Project
from myUtils.TimeTrack import TimeTrack
from myUtils.PlotHelper import PlotHelper
from wingConstruction.wingUtils.defines import *

from mpl_toolkits.mplot3d import axes3d

from matplotlib import cm
from multiprocessing import Pool
import time
from scipy.interpolate import interp1d
from scipy import interpolate

# USED_CORES = 1
# if not USE_ABAQUS:
USED_CORES = Constants().config.getint('meta', 'used_cores', fallback=1)



class MultiRun:

    def __init__(self, use_calcu=True, use_aba=True, non_liner=False, project_name_prefix='pro', force_recalc=False):
        self.use_calculix = use_calcu
        self.use_abaqus = use_aba
        self.non_linear = non_liner
        self.project_name_prefix = project_name_prefix
        self.task_total = 0
        self.task_done = 0
        self.force_recalc = force_recalc

    def print_state(self):
        print('done with {:d} of {:d}'.format(self.task_done, self.task_totals))

    '''
    :param input, a list of the inputs [ribCount, shellThinckness]
    :return only the stress as float
    '''
    def calc_stress(self, input):
        pro = self.new_project_r_t(input[0], input[1])
        pro = self.run_project(pro)
        if self.use_calculix:
            return pro.resultsCalcu.stressMisesMax
        elif self.use_abaqus:
            return pro.resultsAba.stressMisesMax
        return 0.

    def new_project_r_t(self, rib, thick, element_size=0.1):
        if rib % 1 > 0.:
            print('WARNING: rib should be type int but was {:f} in MultiRun.new_project_r_t'.format(rib))
        rib = int(rib)
        project_name = self.project_name_prefix + '_r{:02d}_t{:.6f}'.format(rib, thick)
        pro = self.new_project(project_name)
        pro.ribs = rib
        pro.shellThickness = thick
        pro.elementSize = element_size
        return pro

    def new_project(self, project_name):
        # project_name = 'meshSize_r{:02d}_t{:5f}'.format(rib_count, shell_thick)
        pro = Project(project_name)
        pro.halfSpan = wing_length
        pro.boxDepth = chord_length * 0.4
        pro.boxHeight = chord_height
        pro.ribs = int(wing_length) + 1
        pro.enginePos = engine_pos_y
        pro.engineWeight = engine_weight
        pro.boxOverhang = 0.
        pro.forceTop = -(2. / 3.) * wing_load
        pro.forceBot = -(1. / 3.) * wing_load
        pro.elementSize = 0.1
        # pro1.elementSize = 0.05
        pro.elemType = 'qu4'
        pro.shellThickness = 0.005
        pro.stringerHeight = 0.
        return pro

    def run_project(self, pro):
        if not pro.preexisting or self.force_recalc:
            pro.generate_geometry(nonlinear=self.non_linear)
            if self.use_calculix:
                pro.solve()
                #print('############ DONE ############')
                if not pro.errorFlag:
                    if self.non_linear:
                        pro.post_process(template='wing_post_nl_simple')
                    else:
                        pro.post_process(template='wing_post_simple')
            if self.use_abaqus:
                pro.generate_geometry_abaqus()
                pro.solve_abaqus()
                if not pro.errorFlag:
                    pro.post_process_abaqus()
            pro.save_results()
        #print('#########################################')
        print('finished: ' + pro.workingDir)
        self.task_done += 1
        return pro

    def pool_run(self, projects):
        self.task_done = 0
        self.task_totals = len(projects)
        start = time.time()
        with Pool(USED_CORES) as p:
            projects = p.map(self.run_project, projects)
        print("Time taken = {0:.5f}".format(time.time() - start))
        return projects

    def main_run(self, cleanup=False):
        ribs = np.arange(5, 36, 1)
        ribs = list(ribs)
        thick = np.arange(0.002, 0.0091, 0.0001)
        thick = list(thick)
        projects = []
        for r in ribs:
            for t in thick:
                pro = self.new_project_r_t(r, t)
                projects.append(pro)
        projects = self.pool_run(projects)
        output_file_name = '2drun_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
        output_f = open(Constants().WORKING_DIR + '/'
                        + output_file_name,
                        'w')
        output_f.write(Project.EXPORT_HEADER)
        for p in projects:
            outStr = p.collect_results()
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

    def read_data_file(self, file_name, use_abaqus=False):
        file_path = Constants().WORKING_DIR + '/' + file_name
        data = np.genfromtxt(file_path, delimiter=',', skip_header=1)
        ribs_raw = data[:, 2]
        shell_thick_raw = data[:, 3]
        weight_raw = data[:, 4]
        if not use_abaqus:
            max_stress_raw = data[:, 8]
            max_disp_raw = data[:, 5]
            #max_stress_fixed_raw = data[:, 9]
        else:
            max_stress_raw = data[:, 12]
            max_disp_raw = data[:, 10]

        n_rib = len(set(ribs_raw))
        n_thick = len(set(shell_thick_raw))

        ribs = sorted(list(set(ribs_raw)))
        shell_thick = sorted(list(set(shell_thick_raw)))
        # shellThick = [x * 1000 for x in shellThick]
        weight = np.zeros((n_thick, n_rib))
        max_stress = np.zeros((n_thick, n_rib))
        #max_stress_fixed = np.zeros((n_thick, n_rib))
        max_disp = np.zeros((n_thick, n_rib))
        for i in range(0, len(ribs_raw)):
            weight[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = weight_raw[i]
            max_stress[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = max_stress_raw[i]
            #max_stress_fixed[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = max_stress_fixed_raw[i]
            max_disp[shell_thick.index(shell_thick_raw[i])][ribs.index(ribs_raw[i])] = max_disp_raw[i]

        return ribs, shell_thick, max_stress, max_disp, weight#, max_stress_fixed

    def plot_results(self, file_name):
        import matplotlib.pyplot as plt
        ribs, shell_thick, max_stress, max_disp, weight = self.read_data_file(file_name, use_abaqus=True)
        # max_stress = max_stress_fixed
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
            # if np.min(stress) < max_shear_strength and np.max(stress) > max_shear_strength:
            #   optRibs = np.interp(max_shear_strength, stress, ribs)
            #   f = interp1d(stress, ribs, kind='linear')
            #   plot1.ax.plot([f(max_shear_strength)], [max_shear_strength], 'go')
            #   opti_ribs.append(f(max_shear_strength))
            #   opti_shell.append(shellThick[i])
            plot1.ax.plot(ribs, max_stress[i], label='shell= {:03f}'.format(shell_thick[i]))
        plot1.ax.plot(ribs, np.full((len(ribs), 1), max_shear_strength), 'r--', label='Limit-Load')
        plot1.finalize(legendNcol=2, tighten_layout=False)
        # plot1.show()

        plot2 = PlotHelper(['shellthickness in mm', 'max stress'])
        for i in range(0, len(ribs)):
            stress = max_stress[:, i]
            plot2.ax.plot(shell_thick, stress, label='ribs= {:f}'.format(ribs[i]))
            if np.min(stress) < max_shear_strength and np.max(stress) > max_shear_strength:
                # optRibs = np.interp(max_shear_strength, stress, ribs)
                f = interp1d(stress, shell_thick, kind='linear')
                plot2.ax.plot([f(max_shear_strength)], [max_shear_strength], 'go')
                opti_ribs.append(ribs[i])
                opti_shell.append(f(max_shear_strength))
                opti_weight.append(f_weight(ribs[i], f(max_shear_strength))[0])
        plot2.ax.plot(shell_thick, np.full((len(shell_thick), 1), max_shear_strength), 'r--', label='Limit-Load')
        # plot2.ax.set_xlim((min([x * 1000 for x in shellThick]), max([x * 1000 for x in shellThick])))

        plot2.finalize(legendNcol=2)
        plot2.show()

        opti_min_index = opti_weight.index(min(opti_weight))
        print('opt wight {:03f}'.format(opti_weight[opti_min_index]))
        print('@ {:00f} ribs'.format(opti_ribs[opti_min_index]))
        print('@ {:05f} shell-thickness'.format(opti_shell[opti_min_index]))

        plot3d = PlotHelper(['ribs', 'shell thickness in m', 'mises stress'])
        max_stress[max_stress > 1.2 * max_shear_strength] = np.nan
        color_map = plt.cm.jet(weight / np.max(weight))
        surf = plot3d.ax.plot_surface(rib_mat, shell_mat, max_stress, facecolors=color_map,
                                      cstride=1,
                                      rstride=1)
        optiLine = plot3d.ax.plot(opti_ribs, opti_shell, max_shear_strength, 'k--', label='Limit-Load')
        optiPoint = plot3d.ax.plot([opti_ribs[opti_min_index]], [opti_shell[opti_min_index]], max_shear_strength, 'ro',
                                   label='glob. optimum')
        m = cm.ScalarMappable(cmap=cm.jet)
        m.set_array(weight)
        cbar = plt.colorbar(m)
        cbar.set_label('structure weight in kg', rotation=270)
        limit = np.full((n_thick, n_rib), max_shear_strength)
        plot3d.ax.plot_wireframe(rib_mat, shell_mat, limit, color='r', alpha=0.1)
        plot3d.ax.set_zlim3d(0, max_shear_strength)
        plot3d.finalize()
        plot3d.show()

    def convergence_analysis_run(self, cleanup=False):
        sizes = np.arange(0.08, .26, 0.01)
        sizes = list(sizes)
        projects = []
        for s in sizes:
            pro = self.new_project_r_t(14, shell_thickness, element_size=s)
            projects.append(pro)
        projects = self.pool_run(projects)
        output_file_name = 'convAna_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv'
        output_f = open(Constants().WORKING_DIR + '/' + output_file_name, 'w')
        output_f.write(Project.EXPORT_HEADER)
        for p in projects:
            out_str = p.collect_results()
            if out_str != '':
                output_f.write(out_str)
                output_f.flush()
            else:
                print('ERROR, empty data return')
        output_f.close()
        if cleanup:
            for p in projects:
                p.remove()
        print('DONE with ALL')
        return output_file_name

    def run_sample_points(self, ribs, shells, use_abaqus=False):
        stress = np.zeros((len(ribs)))
        projects = []
        for i in range(0, len(ribs)):
            pro = self.new_project_r_t(int(ribs[i]), shells[i])
            projects.append(pro)
        projects = self.pool_run(projects)
        for i in range(0, len(projects)):
            if projects[i].ribs == ribs[i] and projects[i].shellThickness == shells[i]:
                if use_abaqus:
                    stress[i] = projects[i].resultsAba.stressMisesMax
                else:
                    stress[i] = projects[i].resultsCalcu.stressMisesMax
            else:
                print('ERROR, pool returned shuffled project list')
        return stress


if __name__ == '__main__':
    multi = MultiRun()
    #multi.convergence_analysis_run(cleanup=False)

    # output_file_name = '/dataOut/newRun/2drun_2018-08-10_12_07_48.csv'
    output_file_name = '/dataOut/oldRun/2drun_2018-08-10_12_13_55.csv'
    output_file_name = '2drun_2018-08-22_23_57_54.csv'
    output_file_name = '2drun_2018-08-22_23_05_32.csv'
    output_file_name = '2drun_2018-08-23_16_49_18_final01_cruiseLoad.csv'
    #output_file_name = '2drun_2018-08-23_18_12_26_realLoad25.csv'

    #output_file_name = multi.main_run(cleanup=False)
    multi.plot_results(output_file_name)
