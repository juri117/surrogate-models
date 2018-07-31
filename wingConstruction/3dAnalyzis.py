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

from multiprocessing import Pool
import time

USED_CORES = 4

mtow = 29257
fuel_mass_in_wings = 2*2659.
wing_load = (mtow - fuel_mass_in_wings) * 9.81
engine_mass = 1125 * 9.81
engine_pos_y = 3
wing_length = 12.87
chord_length = 3.
chord_height = 0.55

shear_strength = 3.31e8

element_size = 0.2



def new_project(rib_count, shell_thick):
    projectName = 'meshSize_r{:02d}_t{:5f}'.format(rib_count, shell_thick)
    pro = Project(projectName)
    pro.halfSpan = wing_length
    pro.boxDepth = chord_length * 0.4
    pro.boxHeight = chord_height
    pro.nRibs = rib_count
    pro.boxOverhang = 0.
    pro.forceTop = -0.3 * wing_load
    pro.forceBot = -0.2 * wing_load
    pro.elementSize = element_size
    # pro1.elementSize = 0.05
    pro.elemType = 'qu8'
    pro.shellThickness = shell_thick

    return pro

def run_project(pro):
    pro.generate_geometry(nonlinear=False)
    pro.solve()
    print('############ DONE ############')
    if not pro.errorFlag:
        #pro1.postprocess()
        #if not pro1.errorFlag:
        return True
    return False


def collect_results(pro):
    #projectName = 'meshSize_r{:d}_t{:5f}'.format(rib_count, shell_thick)
    #pro = Project(projectName)
    pro.postprocess(template='wing_post_simple')
    l = pro.validate_load('loadTop.frc')
    l += pro.validate_load('loadBot.frc')
    #l = pro1.validate_load('load.frc')
    loadError = (-0.5*wing_load) - l
    if not pro.errorFlag:
        exportRow = str(element_size) + ',' \
        + str(pro.nRibs) + ',' \
        + str(pro.shellThickness) + ',' \
        + str(pro.clx.dispD3Min) + ','\
        + str(pro.clx.dispD3Max) + ','\
        + str(pro.clx.stressMisesMin) + ','\
        + str(pro.clx.stressMisesMax) + ','\
        + str(loadError)+'\n'
        return exportRow
    return ''


def main_run():
    ribs = np.arange(1, 51, 1)
    ribs = list(ribs)
    thick = np.arange(0.001, 0.021, 0.001)
    thick = list(thick)
    projects = []
    for r in ribs:
        for t in thick:
            projects.append(new_project(r, t))

    start = time.time()
    with Pool(USED_CORES) as p:
        res = p.map(run_project, projects)
    print("Time taken = {0:.5f}".format(time.time() - start))

    outputF = open(Constants().WORKING_DIR + '/'
                   + 'convAna_'
                   + datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
                   + '.csv',
                   'w')
    outputF.write('sizes,nRibs,shellThickness,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,loadError\n')

    for p in projects:
        outStr = collect_results(p)
        if outStr != '':
            outputF.write(outStr)
            outputF.flush()
    outputF.close()
    print("Time taken = {0:.5f}".format(time.time() - start))
    print('DONE with ALL')


if __name__ == '__main__':
    main_run()
