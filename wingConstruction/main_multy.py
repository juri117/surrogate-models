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
max_g = 2.5
safety_fac = 1.5
max_shear_strength = shear_strength * max_g * safety_fac

element_size = 0.2

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
    pro.generate_geometry(nonlinear=False)
    pro.solve()
    print('############ DONE ############')
    if not pro.errorFlag:
        pro.postprocess(template='wing_post_simple')
    return pro

def collect_results(pro):
    l = pro.validate_load('loadTop.frc')
    l += pro.validate_load('loadBot.frc')
    loadError = (-0.5*wing_load) - l
    if not pro.errorFlag:
        exportRow = str(element_size) + ',' \
        + str(pro.nRibs) + ',' \
        + str(pro.shellThickness) + ',' \
        + str(pro.clx.dispD3Min) + ','\
        + str(pro.clx.dispD3Max) + ','\
        + str(pro.clx.stressMisesMin) + ','\
        + str(pro.clx.stressMisesMax) + ',' \
        + str(pro.geo.calc_span_division(pro.halfSpan)) + ',' \
        + str(loadError)+'\n'
        return exportRow
    return ''








def run_test(element_size):
    projectName = 'meshSize_{0:.5f}'.format(element_size)
    pro = Project(projectName)
    pro.halfSpan = 17.
    pro.boxDepth = 2.
    pro.boxHeight = 1.
    pro.ribs = 18
    pro.boxOverhang = 0.1
    pro.forceTop = -0.3*(77000. * 9.81)
    pro.forceBot = -0.2 * (77000. * 9.81)
    #pro.elementSize = 0.25
    pro.elementSize = element_size
    pro.elemType = 'qu8'
    pro.shellThickness = 0.0099
    pro.generate_geometry(nonlinear=True)
    pro.solve()
    print('############ DONE ############')
    if not pro.errorFlag:
        #pro1.postprocess()
        #if not pro1.errorFlag:
        return True
    return False


def collect_results(element_size):
    projectName = 'meshSize_{0:.5f}'.format(element_size)
    pro1 = Project(projectName)
    pro1.postprocess(template='wing_post_simple')
    l = pro1.validate_load('loadTop.frc')
    l += pro1.validate_load('loadBot.frc')
    #l = pro1.validate_load('load.frc')
    loadError = (-0.5 * 77000. * 9.81) - l
    if not pro1.errorFlag:
        exportRow = str(element_size) + ','\
        +str(pro1.clx.dispD3Min) + ','\
        + str(pro1.clx.dispD3Max) + ','\
        + str(pro1.clx.stressMisesMin) + ','\
        + str(pro1.clx.stressMisesMax) + ',' \
        + str(pro1.geo.calc_span_division(pro1.halfSpan)) + ',' \
        + str(loadError)+'\n'
        return exportRow
    return ''


def main_run():
    sizes = np.arange(0.06, .26, 0.01)
    sizes = list(sizes)
    start = time.time()
    with Pool(USED_CORES) as p:
        res = p.map(run_test, sizes)
    print("Time taken = {0:.5f}".format(time.time() - start))

    outputF = open(Constants().WORKING_DIR + '/'
                   + 'convAna_'
                   + datetime.now().strftime('%Y-%m-%d_%H_%M_%S')
                   + '.csv',
                   'w')
    outputF.write('sizes,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,loadError\n')

    for s in sizes:
        outStr = collect_results(s)
        if outStr != '':
            outputF.write(outStr)
            outputF.flush()
    outputF.close()
    print("Time taken = {0:.5f}".format(time.time() - start))
    print('DONE with ALL')


if __name__ == '__main__':
    main_run()
