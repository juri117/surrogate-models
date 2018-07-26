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


t = TimeTrack()
t.tic()
projectName = 'test'
pro1 = Project(projectName)
pro1.halfSpan = 17.
pro1.boxDepth = 2.
pro1.boxHeight = 1.
pro1.nRibs = 30
pro1.boxOverhang = 0.
pro1.forceTop = -0.3*(77000. * 9.81)
pro1.forceBut = -0.2*(77000. * 9.81)
pro1.elementSize = .2
#pro1.elementSize = 0.05
pro1.elemType = 'qu8'
pro1.shellThickness = 0.002
pro1.generate_geometry()
#todo: detect failed mesh generation
pro1.solve()
if not pro1.errorFlag:
    pro1.postprocess()
    if not pro1.errorFlag:
        runTime = t.toc()

        print('min displacement: ' + str(pro1.clx.dispD3Min))
        print('max displacement: ' + str(pro1.clx.dispD3Max))
        print('min mieses stress: ' + str(pro1.clx.stressMisesMin))
        print('max mieses stress: ' + str(pro1.clx.stressMisesMax))

        #l = pro1.validate_load('load.frc')

        l = pro1.validate_load('loadTop.frc')
        l += pro1.validate_load('loadBut.frc')
        print('load error: ' + str((-0.5 * 77000. * 9.81) - l))

