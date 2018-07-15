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

sizes = np.flip(np.arange(0.025, .1, 0.001),0)

outputF = open(Constants().WORKING_DIR + '/' + 'convAna_' + datetime.now().strftime('%Y-%m-%d_%H_%M_%S') + '.csv', 'w')
outputF.write('sizes,dispD3Min,dispD3Max,stressMisesMin,stressMisesMax,time,error\n')


for s in sizes:
    t = TimeTrack()
    t.tic()
    projectName = 'meshSize_{0:.5f}'.format(s)
    pro1 = Project(projectName)
    pro1.halfSpan = 17.
    pro1.boxDepth = 2.
    pro1.boxHeight = 1.
    pro1.nRibs = 17
    pro1.boxOverhang = 0.1
    pro1.forceTop = -0.3*(77000. * 9.81)
    pro1.forceBut = -0.2*(77000. * 9.81)
    #pro1.elementSize = 0.25
    pro1.elementSize = s
    pro1.elemType = 'qu4'
    pro1.shellThickness = 0.0099
    pro1.generate_geometry()
    #todo: detect failed mesh generation
    pro1.solve()
    if not pro1.errorFlag:
        pro1.postprocess()
        if not pro1.errorFlag:
            runTime = t.toc()

            outputF.write(str(s)+','
                          +str(pro1.clx.dispD3Min)+','
                          + str(pro1.clx.dispD3Max)+','
                          + str(pro1.clx.stressMisesMin)+','
                          + str(pro1.clx.stressMisesMax)+','
                          + str(runTime)+','
                          + str(pro1.errorFlag)+'\n')
            outputF.flush()
            print('min displacement: ' + str(pro1.clx.dispD3Min))
            print('max displacement: ' + str(pro1.clx.dispD3Max))
            print('min mieses stress: ' + str(pro1.clx.stressMisesMin))
            print('max mieses stress: ' + str(pro1.clx.stressMisesMax))

outputF.close()
