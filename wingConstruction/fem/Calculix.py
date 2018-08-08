__author__ = "Juri Bieler"
__version__ = "0.0.1"
__status__ = "Development"

# ==============================================================================
# description     :interface to Calculix and solver-input-file generation
# author          :Juri Bieler
# date            :2018-07-13
# notes           :
# python_version  :3.6
# ==============================================================================

import subprocess
import os

from wingConstruction.utils.Constants import Constants

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
        self._workingDir = workingDir
        if not os.path.isdir(self._workingDir):
            os.mkdir(self._workingDir)
        self.dispD3Min = 0
        self.dispD3Max = 0
        self.stressMisesMin = 0
        self.stressMisesMax = 0
        self.stressMisesMaxFixed = 0

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
    def generate_mesh(self, part_name):
        print('run fem pre-processor cgx('+self._workingDir+')')
        p = self.run_cgx(part_name + '.fbl', pipe_response=True)
        out, err = p.communicate()

    # calls the fem solver (all input files must be present in the working directory)
    def solve_model(self, file_name):
        # print('--- start ccx output ---------------------------------------')
        print('run fem solver ccx('+self._workingDir+')')
        p = self.run_ccx(file_name, pipe_response=True)
        out, err = p.communicate()
        """
        out, err = p.communicate()
        print(out.decode('UTF-8'))
        # print('--- end ccx output ---------------------------------------')
        if err is not None and err != b'':
            self._log.error('ccx process failed')
            self.errorFlag = True
        if 'ERROR' in out.decode('UTF-8'):
            print('ccx returned an error')
            self.errorFlag = True
        if 'Job finished' in out.decode('UTF-8'):
            print('ccx succeeded')
        """

    def run_postprocessing(self, file_name):
        # ToDo: run it with pipeRespnose=True results in calculix gui freeze
        # print('--- start cgx output ---------------------------------------')
        print('run fem post-processing cgx('+self._workingDir+')')
        p = self.run_cgx(file_name, pipe_response=True)
        out, err = p.communicate()
        #print(out.decode('UTF-8'))
        # print('--- stop cgx output ---------------------------------------')
        if err is not None and err != b'':
            print('ccx process failed('+self._workingDir+')')
            self.errorFlag = True
        elif 'ERROR' in out.decode('UTF-8'):
            print('ccx returned an error('+self._workingDir+')')
            self.errorFlag = True
        else:
            print('cgx succeeded('+self._workingDir+')')
            self.process_cgx_output(out)

    def process_cgx_output(self, strOut):
        lines = strOut.split(b'\n')
        for i in range(0, len(lines) - 3):
            if b'DISP' in lines[i] and b'D3' in lines[i]:
                if b'max:' in lines[i + 2]:
                    str = lines[i + 2].replace(b' max:', b'')
                    str = str.split(b' ')[0]
                    self.dispD3Max = float(str)
                if b'min:' in lines[i + 3]:
                    str = lines[i + 3].replace(b' min:', b'')
                    str = str.split(b' ')[0]
                    self.dispD3Min = float(str)
            if b'STRESS' in lines[i] and b'Mises' in lines[i]:
                if b'max:' in lines[i + 2]:
                    str = lines[i + 2].replace(b' max:', b'')
                    str = str.split(b' ')[0]
                    self.stressMisesMax = float(str)
                if b'min:' in lines[i + 3]:
                    str = lines[i + 3].replace(b' min:', b'')
                    str = str.split(b' ')[0]
                    self.stressMisesMin = float(str)
            if b'mode:h' in lines[i]:
                if b'node' in lines[i+1] and b'value' in lines[i+1]:
                    str = lines[i + 1].replace(b' node:', b'')
                    str = str.replace(b'value:', b'')
                    str = str.replace(b'dist:', b'')
                    vals = str.split(b' ')
                    self.stressMisesMaxFixed = float(vals[1])

    ##############################################
    # helper functions

    def run_ccx(self, fileName, pipe_response=False):
        if pipe_response:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self._workingDir,
                                 stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CCX_EXE_PATH, fileName], cwd=self._workingDir)
            p.wait()
        return p

    def run_cgx(self, file_name, pipe_response=False):
        if pipe_response:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', file_name], cwd=self._workingDir,
                                 stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        else:
            p = subprocess.Popen([Constants().CALCULIX_CGX_EXE_PATH, '-b', file_name], cwd=self._workingDir)
            p.wait()
        return p

    '''
    def run_cgx_postprocessing(self, partName):
        p = subprocess.Popen([self.pathCGX, '-b', partName + '.fbd'], cwd=self.workingDir)
        p.wait()
        print("done")
        return p
    '''

    def get_file_path(self, fileName):
        return self._workingDir + '/' + fileName


if __name__ == '__main__':
    clx = Calculix(workingDir='../dataOut/test01')
    clx.generate_mesh('test')
    #clx.generate_geometry()
    #clx.generate_mesh()
    #clx.generate_bc()
    #clx.generate_load()
    #clx.generate_inp()
    #clx.solve_model()
