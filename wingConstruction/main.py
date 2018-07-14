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

from wingConstruction.Project import Project

pro1 = Project('firstTry')
pro1.halfSpan = 17.
pro1.boxDepth = 2.
pro1.boxHeight = 1.
pro1.nRibs = 9
pro1.boxOverhang = 0.1
pro1.forceTop = -0.3*(77000. * 9.81)
pro1.forceBut = -0.2*(77000. * 9.81)
pro1.elementSize = 0.25
pro1.elemType = 'qu8'
pro1.shellThickness = 0.0099
pro1.generate_geometry()
pro1.solve()
pro1.postprocess()

