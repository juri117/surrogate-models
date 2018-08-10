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


max_g = 2.5
safety_fac = 1.5
mtow = 29257
fuel_mass_in_wings = 2*2659.
wing_load = ((mtow - fuel_mass_in_wings - 1000) * 9.81) * max_g * safety_fac
engine_weight = 1125. * 9.81
engine_pos_y = 3.
wing_length = 12.87
chord_length = 3.
chord_height = 0.55

shear_strength = 3.31e8

t = TimeTrack()
t.tic()
projectName = 'test_14_eng'
pro1 = Project(projectName)
pro1.halfSpan = wing_length
pro1.boxDepth = chord_length*0.4
pro1.boxHeight = chord_height
pro1.ribs = 14
pro1.enginePos = engine_pos_y
pro1.engineWeight = engine_weight
pro1.boxOverhang = 0.
pro1.forceTop = -0.3*wing_load
pro1.forceBot = -0.2*wing_load
pro1.elementSize = .1
#pro1.elementSize = 0.05
pro1.elemType = 'qu4'
pro1.shellThickness = 0.009
pro1.generate_geometry(nonlinear=False)
#todo: detect failed mesh generation
pro1.solve()
if not pro1.errorFlag:
    pro1.post_process()
    if not pro1.errorFlag:
        runTime = t.toc()

        print('min displacement: ' + str(pro1.clx.dispD3Min))
        print('max displacement: ' + str(pro1.clx.dispD3Max))
        print('min mieses stress: ' + str(pro1.clx.stressMisesMin))
        print('max mieses stress: ' + str(pro1.clx.stressMisesMax))

        #l = pro1.validate_load('load.frc')

        l = pro1.validate_load('loadTop.frc')
        l += pro1.validate_load('loadBot.frc')
        print('load error: ' + str((-0.5*wing_load) - l))

