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
import sys
from datetime import datetime
from wingConstruction.utils.Constants import Constants
from wingConstruction.Project import Project
from utils.TimeTrack import TimeTrack
from wingConstruction.utils.defines import *

t = TimeTrack()
t.tic()
projectName = 'test000'#_r07_s002'
pro1 = Project(projectName)
pro1.halfSpan = wing_length
pro1.boxDepth = chord_length*0.4
pro1.boxHeight = chord_height
pro1.ribs = 18
pro1.enginePos = engine_pos_y
pro1.engineWeight = engine_weight
pro1.boxOverhang = 0.
pro1.forceTop = -(2./3.) * wing_load
pro1.forceBot = -(1./3.) * wing_load
pro1.elementSize = .1
#pro1.elementSize = 0.05
pro1.elemType = 'qu4'
pro1.shellThickness = 0.008
pro1.stringerHeight = 0.
pro1.generate_geometry(nonlinear=False)

pro1.generate_geometry_abaqus()
pro1.solve_abaqus()
pro1.post_process_abaqus()

#todo: detect failed mesh generation
pro1.solve()
if not pro1.errorFlag:
    pro1.post_process()
    #pro1.post_process(template='wing_post_max_mises_fixed')
    if not pro1.errorFlag:
        runTime = t.toc()

        print('min displacement: ' + str(pro1.clx.dispD3Min))
        print('max displacement: ' + str(pro1.clx.dispD3Max))
        print('min mieses stress: ' + str(pro1.clx.stressMisesMin))
        print('max mieses stress: ' + str(pro1.clx.stressMisesMax))
        print('fix mieses stress: ' + str(pro1.clx.stressMisesMaxFixed))

        #l = pro1.validate_load('load.frc')

        l = pro1.validate_load('loadTop.frc')
        l += pro1.validate_load('loadBot.frc')
        print('load error: ' + str((-1.*wing_load) - l))

print('done')