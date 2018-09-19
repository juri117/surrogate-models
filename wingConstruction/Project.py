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
import numpy as np

from wingConstruction.utils.Constants import Constants
from wingConstruction.fem.WingConstructionV4 import WingConstruction
from wingConstruction.fem.Calculix import Calculix
from wingConstruction.fem.Abaqus import Abaqus


class Project:

    EXPORT_HEADER = 'elementSizes,spanElementCount,ribs,shellThickness,weight,' \
                    + 'dispD3Min(Cal),dispD3Max(Cal),stressMisesMin(Cal),stressMisesMax(Cal),' \
                    + 'dispD3Min(Aba),dispD3Max(Aba),stressMisesMin(Aba),stressMisesMax(Aba),loadError\n'

    def __init__(self, project_name):
        self.errorFlag = False
        # indicates if the project was calculated before
        self.preexisting = False

        self.resultsCalcu = ResultMax()
        self.resultsAba = ResultMax()

        self.workingDir = Constants().WORKING_DIR + '/' + project_name
        if not os.path.isdir(self.workingDir):
            os.mkdir(self.workingDir)
        else:
            if os.path.isfile(self.workingDir + '/' + 'results.csv'):
                self.preexisting = True
                self.parse_from_results()
                # hack to recalculate failed projects
                if self.resultsCalcu.stressMisesMax == 0 and self.resultsAba.stressMisesMax == 0:
                    self.preexisting = False

        self.clx = Calculix(workingDir=self.workingDir)
        self.geo = None
        self.aba = None

        self.halfSpan = 17.
        self.boxDepth = 2.
        self.boxHeight = 1.
        self.ribs = 8
        self.enginePos = 3.
        self.engineWeight = 1000. * 9.81
        self.boxOverhang = 0.5
        self.stringerHeight = 0.
        self.forceTop = 0.
        self.forceBot = 0.
        self.elementSize = 0.25
        self.elemType = 'qu4'
        self.shellThickness = 0.01
        self.density = 2810  # kg/m^3

    def _get_geo(self):
        if self.geo is None:
            self.geo = WingConstruction(self.workingDir,
                                        self.halfSpan,
                                        self.boxDepth,
                                        self.boxHeight,
                                        self.ribs,
                                        self.shellThickness,
                                        self.enginePos,
                                        box_overhang=self.boxOverhang,
                                        stringer_height=self.stringerHeight)
            self.geo.element_size = self.elementSize
        return self.geo

    def generate_geometry(self, nonlinear=False):
        self._get_geo().generate_wing(self.forceTop,
                               self.forceBot,
                               self.engineWeight,
                               self.elementSize,
                               element_type=self.elemType)
        self._get_geo().generate_inp(nonlinear=nonlinear)
        if self.clx is None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.generate_mesh('wingGeo')

    def calc_span_division(self):
        return self._get_geo().calc_span_division(self.halfSpan)

    def generate_geometry_abaqus(self):
        if self.aba is None:
            self.aba = Abaqus(self.workingDir)
        self.aba.calculix_to_abaqus(self.shellThickness)

    def solve_abaqus(self):
        if self.aba is None:
            self.aba = Abaqus(self.workingDir)
        self.aba.solve_model()
        if self.aba.errorFlag:
            self.errorFlag = True

    def post_process_abaqus(self):
        if self.aba is None:
            self.aba = Abaqus(self.workingDir)
        self.aba.post_processing()
        self.resultsAba.dispD3Min = self.aba.dispD3Min
        self.resultsAba.dispD3Max = self.aba.dispD3Max
        self.resultsAba.stressMisesMin = self.aba.stressMisesMin
        self.resultsAba.stressMisesMax = self.aba.stressMisesMax

    def solve(self):
        if self.clx is None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.solve_model('wingRun')
        if self.clx.errorFlag:
            self.errorFlag = True

    def post_process(self, template='wing_post'):
        copyfile(Constants().INPUT_DIR + '/' + template + '.fbd', self.workingDir + '/'+template+'.fbd')
        if self.clx is None:
            self.clx = Calculix(workingDir=self.workingDir)
        self.clx.run_postprocessing(template+'.fbd')
        if self.clx.errorFlag:
            self.errorFlag = True
        self.resultsCalcu.dispD3Min = self.clx.dispD3Min
        self.resultsCalcu.dispD3Max = self.clx.dispD3Max
        self.resultsCalcu.stressMisesMin = self.clx.stressMisesMin
        self.resultsCalcu.stressMisesMax = self.clx.stressMisesMax

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

    def calc_wight(self):
        return self._get_geo().calc_weight(self.density)

    def collect_results(self):
        l = self.validate_load('loadTop.frc')
        l += self.validate_load('loadBot.frc')
        load_error = (self.forceTop + self.forceBot) - l
        export_row = str(self.elementSize) + ',' \
                     + str(self.calc_span_division()) + ',' \
                     + str(self.ribs) + ',' \
                     + str(self.shellThickness) + ',' \
                     + str(self.calc_wight()) + ',' \
                     + str(self.resultsCalcu.dispD3Min) + ',' \
                     + str(self.resultsCalcu.dispD3Max) + ',' \
                     + str(self.resultsCalcu.stressMisesMin) + ',' \
                     + str(self.resultsCalcu.stressMisesMax) + ',' \
                     + str(self.resultsAba.dispD3Min) + ',' \
                     + str(self.resultsAba.dispD3Max) + ',' \
                     + str(self.resultsAba.stressMisesMin) + ',' \
                     + str(self.resultsAba.stressMisesMax) + ',' \
                     + str(load_error) + '\n'
        return export_row

    def save_results(self):
        output_f = open(self.workingDir + '/' + 'results.csv', 'w')
        output_f.write(Project.EXPORT_HEADER)
        output_f.write(self.collect_results())
        output_f.close()

    def parse_from_results(self):
        data = np.genfromtxt(self.workingDir + '/' + 'results.csv', delimiter=',', skip_header=1)
        self.resultsCalcu.stressMisesMin = data[7]
        self.resultsCalcu.stressMisesMax = data[8]
        self.resultsCalcu.dispD3Min = data[6]
        self.resultsCalcu.dispD3Max = data[5]
        self.resultsAba.stressMisesMin = data[11]
        self.resultsAba.stressMisesMax = data[12]
        self.resultsAba.dispD3Min = data[9]
        self.resultsAba.dispD3Max = data[10]

    def remove(self):
        rmtree(self.workingDir)


class ResultMax:
    def __init__(self):
        self.dispD3Min = 0
        self.dispD3Max = 0
        self.stressMisesMin = 0
        self.stressMisesMax = 0
