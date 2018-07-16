__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :creates a new project and offers all functionality to run it
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import os
from shutil import copyfile
from shutil import rmtree

from wingConstruction.utils.Constants import Constants
from wingConstruction.fem.WingConstruction import WingConstruction
from wingConstruction.fem.Calculix import Calculix

class Project:

    def __init__(self, project_name):
        self.errorFlag = False
        self.workingDir = Constants().WORKING_DIR+'/'+project_name
        if not os.path.isdir(self.workingDir):
            os.mkdir(self.workingDir)

        self.clx = None

        self.halfSpan = 17.
        self.boxDepth = 2.
        self.boxHeight = 1.
        self.nRibs = 8
        self.boxOverhang = 0.5
        self.forceTop = 0.
        self.forceBut = 0.
        self.elementSize = 0.25
        self.elemType = 'qu4'
        self.shellThickness = 0.01


    def generate_geometry(self):
        geo = WingConstruction(self.workingDir,
                               self.halfSpan,
                               self.boxDepth,
                               self.boxHeight,
                               self.nRibs,
                               box_overhang=self.boxOverhang)

        geo.generate_wing(self.forceTop,
                          self.forceBut,
                          self.elementSize,
                          element_type=self.elemType)
        geo.generate_inp(self.shellThickness)

    def solve(self):
        if self.clx == None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.generate_mesh('wingGeo')
        self.clx.solve_model('wingRun')
        if self.clx.errorFlag:
            self.errorFlag = True

    def postprocess(self, template='wing_post'):
        copyfile(Constants().INPUT_DIR+'/'+template+'.fbd', self.workingDir+'/wing_post.fbd')
        if self.clx == None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.run_postprocessing('wing_post.fbd')
        if self.clx.errorFlag:
            self.errorFlag = True

    def remove(self):
        rmtree(self.workingDir)