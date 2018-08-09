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
from wingConstruction.fem.WingConstructionV3 import WingConstruction
from wingConstruction.fem.Calculix import Calculix


class Project:

    def __init__(self, project_name):
        self.errorFlag = False
        self.workingDir = Constants().WORKING_DIR + '/' + project_name
        if not os.path.isdir(self.workingDir):
            os.mkdir(self.workingDir)

        self.clx = None
        self.geo = None

        self.halfSpan = 17.
        self.boxDepth = 2.
        self.boxHeight = 1.
        self.ribs = 8
        self.boxOverhang = 0.5
        self.forceTop = 0.
        self.forceBot = 0.
        self.elementSize = 0.25
        self.elemType = 'qu4'
        self.shellThickness = 0.01

    def generate_geometry(self, nonlinear=False):
        if self.geo is None:
            self.geo = WingConstruction(self.workingDir,
                                        self.halfSpan,
                                        self.boxDepth,
                                        self.boxHeight,
                                        self.ribs,
                                        self.shellThickness,
                                        box_overhang=self.boxOverhang)
        self.geo.generate_wing(self.forceTop,
                               self.forceBot,
                               self.elementSize,
                               element_type=self.elemType)
        self.geo.generate_inp(nonlinear=nonlinear)

    def solve(self):
        if self.clx is None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.generate_mesh('wingGeo')
        self.clx.solve_model('wingRun')
        if self.clx.errorFlag:
            self.errorFlag = True

    def post_process(self, template='wing_post'):
        copyfile(Constants().INPUT_DIR + '/' + template + '.fbd', self.workingDir + '/wing_post.fbd')
        if self.clx is None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.run_postprocessing('wing_post.fbd')
        if self.clx.errorFlag:
            self.errorFlag = True

    def validate_load(self, load_file_name):
        load_f = open(self.workingDir + '/' + load_file_name)
        load_sum = 0
        for line in load_f:
            try:
                values = line.strip().split(',')
                if len(values) == 3:
                    load = float(values[2].strip())
                    load_sum += load
            except ValueError:
                print('invalid line in Load file')
        print('sum of Loads in ' + load_file_name + ': ' + str(load_sum))
        return load_sum

    def remove(self):
        rmtree(self.workingDir)
