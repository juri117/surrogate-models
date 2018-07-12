# ==============================================================================
# description     :interface to Calculix and solver-input-file generation
# author          :Juri Bieler
# date            :2017-12-13
# version         :0.01
# notes           :
# python_version  :3.6
# ==============================================================================

# plot:
# plot la all -> plot all nodes with index as label

import subprocess
import os
import matplotlib.pyplot as plt
import traceback
import math
from time import sleep

from Constants import Constants

##############################################
# set up variables

#SHELL_THICKNESS = 0.008 #in mm
#SHELL_MATERIAL_YOUNG = 60e9 #N/m^2
#SHELL_MATERIAL_POISSON = 0.34

class Calculix():

    ##############################################
    # system variables [DO NOT TOUCH]
    pointCount = 0
    lineCount = 0
    bcs = []
    errorFlag = False

    # init function prepares working environment
    def __init__(self, workingDir='data'):
        self.workingDir = workingDir
        if not os.path.isdir(self.workingDir):
            os.mkdir(self.workingDir)

    ##############################################
    # fem postprocessing

    def generate_postpro_file(self, partName):
        f = open(self.get_file_path(partName + '_post.fbd'), 'w')
        inpLines = []
        inpLines.append('read ' + partName + '.frd')
        inpLines.append('read ' + partName + '.inp nom') #nom -> no-mesh, since it is already in frd file
        # deformation
        inpLines.append('view edge off')
        inpLines.append('text deformation')
        inpLines.append('rot y')
        inpLines.append('rot r 90')
        inpLines.append('ds 1 e 3')
        inpLines.append('view disp')
        #inpLines.append('seta ! all')
        inpLines.append('frame')
        inpLines.append('hcpy png')
        inpLines.append('sys mv hcpy_1.png ' + partName + '_deform.png')
        # mieses top
        inpLines.append('text Mieses stress top')
        inpLines.append('view disp off')
        inpLines.append('rot y')
        inpLines.append('rot u 90')
        inpLines.append('ds 2 e 7')
        inpLines.append('min 0 f')
        inpLines.append('frame')
        inpLines.append('hcpy png')
        inpLines.append('sys mv hcpy_2.png ' + partName + '_mieses_top.png')
        #mieses buttom
        inpLines.append('text Mieses stress buttom')
        inpLines.append('rot r 180')
        inpLines.append('frame')
        inpLines.append('hcpy png')
        inpLines.append('sys mv hcpy_3.png ' + partName + '_mieses_buttom.png')
        #mesh buttom
        inpLines.append('view surf')
        inpLines.append('view elem')
        inpLines.append('frame')
        inpLines.append('hcpy png')
        inpLines.append('sys mv hcpy_4.png ' + partName + '_mesh_buttom.png')
        inpLines.append('quit')

        f.writelines(line + '\n' for line in inpLines)
        f.close()

    ##############################################
    # fem calls

    # generates the mesh from a given .fbl-file
    def generate_mesh(self, partName):
        p = self.run_cgx(partName + '.fbl')

    # calls the fem solver (all input files must be present in the working directory)
    def solve_model(self, partName):
        p = self.run_ccx(partName, pipeResponse=False)

    def run_postprocessing(self, partName):
        #ToDo: run it with pipeRespnose=True results in calculix gui freeze
        print('--- start cgx output ---------------------------------------')
        p = self.run_cgx(partName + '_post.fbd', pipeResponse=False)
        print('--- stop cgx output ---------------------------------------')
        #out, err = p.communicate()
        #print(out.decode('UTF-8'))

    ##############################################
    # helper functions

    def run_ccx(self, fileName, pipeResponse=False):
        if pipeResponse:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self.workingDir, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self.workingDir)
        p.wait()
        print("done")
        return p

    def run_cgx(self, fileName, pipeResponse=False):
        if pipeResponse:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', fileName], cwd=self.workingDir, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', fileName], cwd=self.workingDir)
        p.wait()
        print("done")
        return p

    '''
    def run_cgx_postprocessing(self, partName):
        p = subprocess.Popen([self.pathCGX, '-b', partName + '.fbd'], cwd=self.workingDir)
        p.wait()
        print("done")
        return p
    '''

    def get_file_path(self, fileName):
        return self.workingDir + '/' + fileName


if __name__ == '__main__':
    clx = Calculix(workingDir='../dataOut/test01')
    clx.generate_mesh('test')
    #clx.generate_geometry()
    #clx.generate_mesh()
    #clx.generate_bc()
    #clx.generate_load()
    #clx.generate_inp()
    #clx.solve_model()
